# switchpattern defaults to first 4 days 0-30, 2 days 30 to 70, last 5 days 70+
#' @import mgcv
#' @importFrom utils head
#' @keywords internal
#' @export
prepCumulativeEffectDiffs <- function(
  ped,
  m,
  id,
  which.patients = c("low", "mid", "full", "switch"),
  effectname     = "AdequacyCals",
  switchpattern  = c(rep("low", 4), rep("mid", 2), rep("full", 5))) {

  gc <- grep("AdequacyCals", names(coef(m)), value = TRUE)

  if (length(gc) == 0) {
    effectname <- "ProtCat"

    low.var    <- paste0(effectname, "1")
    mid.var    <- paste0(effectname, "2")
    full.var   <- paste0(effectname, "3")
  } else {
    effectname <- sub(".*I\\(", "", gc)
    effectname <- unique(sub("Tot.*", "", effectname))

    low.var    <- paste0(effectname, "Tot0to30")
    mid.var    <- paste0(effectname, "Tot30To70")
    full.var   <- paste0(effectname, "TotAbove70")
  }

  # pick a patient that remained under risk all the way to overwrite their data:
  protoPat <- {
    ind <- which(ped$intmid > 58)[1]
    id <- ped$CombinedID[ind]
    subset(ped, CombinedID == id)
  }
  pat0 <- patMid <- patFull <- patSwitch <- protoPat

  maxdays.nutri <- ncol(ped$DaysMat)

  stopifnot(length(switchpattern) == maxdays.nutri,
              all(switchpattern %in% c("low", "mid", "full")))

  ## load model if m character
  if( is.character(m) ) m <- readRDS(m)

  ##
  mat.full  <- matrix(1, nrow=nrow(protoPat), ncol=maxdays.nutri)
  mat.empty <- matrix(0, nrow=nrow(protoPat), ncol=maxdays.nutri)


  ## low
  pat0[[low.var]] <- mat.full
  pat0[[mid.var]] <- pat0[[full.var]] <- mat.empty
  ## mid
  patMid[[low.var]]  <- mat.empty
  patMid[[mid.var]]  <- mat.full
  patMid[[full.var]] <- mat.empty
  ## full
  patFull[[low.var]]  <- mat.empty
  patFull[[mid.var]]  <- mat.empty
  patFull[[full.var]] <- mat.full
  # switch
  patSwitch[[low.var]]  <- outer(rep(1, nrow(protoPat)), switchpattern == "low")
  patSwitch[[mid.var]]  <- outer(rep(1, nrow(protoPat)), switchpattern == "mid")
  patSwitch[[full.var]] <- outer(rep(1, nrow(protoPat)), switchpattern == "full")

  diffs <- expand.grid(
    patient1=ordered(which.patients, c("low", "mid", "switch", "full")),
    patient2=ordered(which.patients, c("low", "mid", "switch", "full")))
  diffs <- subset(diffs, patient1 < patient2)

  switchProtoPat <- function(pat){
    switch(as.character(pat),
           "low"    = pat0,
           "mid"    = patMid,
           "switch" = patSwitch,
           "full"   = patFull)
  }
  pat1 <- lapply(diffs[, 1], switchProtoPat)
  pat2 <- lapply(diffs[, 2], switchProtoPat)

  # apply getEffectDiffs over diffs
  ## mapply stupidly copies m around for this, takes 4eva.
  #    diff.list <- mapply(getEffectDiffs, newdata1=pat1, newdata2=pat2,
  #                        MoreArgs=list(m=m, effectname=effectname))
  diff.list <- list()
  for(i in seq_len(nrow(diffs))){

    suppressWarnings(
      diff.list <- append(diff.list, list(cbind(
        diffs[i,],
        getEffectDiffs(
          newdata1   = pat1[[i]],
          newdata2   = pat2[[i]],
          m          = m,
          effectname = effectname)))))

  }
  diff <- Reduce(rbind, diff.list)

  switchpattern_string <- {
    rles      <- rle(switchpattern)
    start     <- ifelse(maxdays.nutri==12, 0, 1)
    startdays <- head(c(start, start + cumsum(rles$lengths)), -1)
    stopdays  <- start + cumsum(rles$lengths) - 1
    tmp       <- paste0(startdays, "-", stopdays, ": ", rles$values)
    tmp[1]    <- paste("days", tmp[1])
    lvls      <- if (grepl("AdequacyCals", effectname)) {
      c("C I", "C II ", "C III")
    } else {
      c("<0.6 g/kg", "0.6-1.2 g/kg", ">1.2 g/kg")
    }
    paste(
        gsub(pattern="low",  replacement=lvls[1],
        gsub(pattern="mid",  replacement=lvls[2],
        gsub(pattern="full", replacement=lvls[3], x=tmp))),
      collapse="; ")
  }


  if (grepl("AdequacyCals", effectname)) {
    patterns <- c(
      "days 1-11: C I ",
      "days 1-11: C II",
      switchpattern_string,
      "days 1-11: C III")
  } else {
    patterns <- c(
      "days 1-11: < 0.6 g/kg",
      "days 1-11: 0.6 - 1.2 g/kg",
      switchpattern_string,
      "days 1-11: > 1.2 g/kg")
  }
  diff$patient1_simple <- diff$patient1
  diff$patient2_simple <- diff$patient2
  levels(diff$patient1) <- levels(diff$patient2) <- patterns

  diff

}


#' Extract effect differences for 2 hypothetical patients
#'
#' @param effectname Name of variable(s) for which the effect should be extracted.
#' @import mgcv
#' @importFrom stats predict
#' @keywords internal
#' @export
getEffectDiffs <- function(
  m,
  newdata1,
  newdata2,
  effectname = "AdequacyCals") {

  ## deal with fits created with old mgcv before ti() was implemented
  ## when loaded with newer versions that have ti(): add the missing "mc"-slot
  if(!is.null(mgcv::ti)){
    m$smooth <- lapply(m$smooth, function(trm){
      if("margin" %in% names(trm) & is.null(trm$mc)){
        trm$mc <- rep(TRUE, length(trm$margin))
                #message("added mc slot to ", trm$label)
      }
      trm
    } )
  }

  useColumns <- grep(effectname, names(m$coefficients))
  covCoefs <- m$Vp[useColumns, useColumns]
  X1 <- predict(m, newdata = newdata1, type = "lpmatrix")[,useColumns]
  X2 <- predict(m, newdata = newdata2, type = "lpmatrix")[,useColumns]
  X <- X2 - X1
  fit <- drop(X %*% m$coefficients[useColumns])
  se <- sqrt(rowSums((X %*% covCoefs) * X))
  hi <- fit + 2 * se
  lo <- fit - 2 * se
  cbind(intmid = newdata1$intmid, fit = fit, se = se, hi = hi, lo = lo)

}


#' @import ggplot2
#' @keywords internal
#' @export
plotCumulativeNutritionDiffs <- function(
  diffs,
  ncol = 2,
  nrow =  ceiling(nlevels(diffs$label)/ncol)) {
  ## interpretation: vergleich der hazardraten zwischen erster
  ## ernaehrung und zweiter ernaehrung
  ## wobei der geplottete effekt "zweite - erste" ist,
  ## also: positive werte bedeuten mehr Risiko fuer "zweite",
  ##       negative Werte bedeuten mehr Risiko fuer "erste"

  p <- ggplot(diffs) +
    geom_hline(yintercept = 0) +
    geom_pointrange(aes(x = intmid, y = fit, ymin = lo, ymax = hi, col = Model),
      #group = as.factor(patient2), color = patient2),
      position = position_dodge(width = .8)) +
    ylab("log-hazard difference") +
    facet_wrap(~label, nrow=2) +
    xlab("Days after ICU admission") +
    theme(legend.position = "bottom")

  p

}
