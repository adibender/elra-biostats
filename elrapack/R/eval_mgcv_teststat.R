#' Label patterns
#'
#' @importFrom utils head
#' @keywords internal
#' @export
pattern_label <- function(
    pattern,
    maxdays.tdc,
    effectname="AdequacyCals") {

    patt.string <- {
        rles      <- rle(as.character(pattern))
        start     <- ifelse(maxdays.tdc==12, 0, 1)
        startdays <- head(c(start, start + cumsum(rles$lengths)), -1)
        stopdays  <- start + cumsum(rles$lengths) - 1
        tmp       <- paste0(startdays, "-", stopdays, ": ", rles$values)
        tmp[1]    <- paste("days", tmp[1])
        lvls      <- if(grepl("AdequacyCals", effectname)) {
            c("C I", "C II", "C III")
        } else {
            c("<0.6 g/kg", "0.6-1.2 g/kg", ">1.2 g/kg")
        }
        paste(
            gsub(pattern="low",  replacement=lvls[1],
                gsub(pattern="mid",  replacement=lvls[2],
                    gsub(pattern="full", replacement=lvls[3], x=tmp))),
            collapse="; ")
    }

    # return
    patt.string
}



#' get differences of effects on log hazard
#'
#' @importFrom reshape2 melt
#' @keywords internal
#' @export
get_comp_diffs <- function(
  patients.list,
  protocols.df,
  protocols.to.compare,
  model,
  effectname = "AdequacyCals") {

  maxdays.tdc <- nrow(protocols.df)

  diff.labs <- apply(protocols.to.compare, 1, function(z) {
    paste(pattern_label(protocols.df[, z[1]], maxdays.tdc=maxdays.tdc), "vs. \n",
      pattern_label(protocols.df[, z[2]], maxdays.tdc=maxdays.tdc))
  })

  effect.diffs <- apply(protocols.to.compare, 1, function(z) {
    diff.i <- as.data.frame(
      getEffectDiffs(
        m=model,
        patients.list[[z[1]]],
        patients.list[[z[2]]],
        effectname = effectname))
  })

  names(effect.diffs) <- diff.labs
  m.effect.diffs <- melt(effect.diffs, id.vars=c("fit", "intmid", "se", "hi", "lo"))
  m.effect.diffs$label <- factor(m.effect.diffs$L1, levels = diff.labs,
    labels=diff.labs)
  m.effect.diffs$comparison <- factor(m.effect.diffs$L1,
    levels=diff.labs, labels=levels(protocols.to.compare[, "comparison"]))

  ## return
  m.effect.diffs

}

#' @keywords internal
#' @export
pattern_pat <- function(
  pattern,
  ped,
  m,
  median.patient,
  effectname    = "AdequacyCals") {

  gc <- grep(effectname, names(coef(m)), value = TRUE)

  effectname <- sub(".*I\\(", "", gc)
  effectname <- unique(sub("Tot.*", "", effectname))

  low.var    <- paste0(effectname, "Tot0to30")
  mid.var    <- paste0(effectname, "Tot30To70")
  full.var   <- paste0(effectname, "TotAbove70")

  # pick a patient that remained under risk all the way to overwrite their data:
  protoPat <- {
    ind <- which(ped$intmid == max(ped$intmid))[1]
    id <- ped$CombinedID[ind]
    subset(ped, CombinedID == id)
  }
  # overwrite with median data for confounders:
  for(var in colnames(median.patient)){
    protoPat[, var] <- median.patient[1, var]
  }

  maxdays.tdc <- ncol(ped$DaysMat)
  max_per_int <- rowSums(protoPat$LHartlDynf)

  stopifnot(length(pattern) == maxdays.tdc,
    all(pattern %in% c("low", "mid", "full")))

  # switch
  if(is.matrix(protoPat[[low.var]])) {
    protoPat[[low.var]]  <- outer(rep(1, nrow(protoPat)), pattern == "low")
    protoPat[[mid.var]]  <- outer(rep(1, nrow(protoPat)), pattern == "mid")
    protoPat[[full.var]] <- outer(rep(1, nrow(protoPat)), pattern == "full")
  } else {
    if(effectname == "LLadequ") {
        protoPat[[low.var]] <- c(0, (rowSums(outer(rep(1, nrow(protoPat)),
          pattern == "low")*protoPat$LHartlDynf)/max_per_int)[-1])
        protoPat[[mid.var]] <- c(0,(rowSums(outer(rep(1, nrow(protoPat)),
          pattern == "mid")*protoPat$LHartlDynf)/max_per_int)[-1])
        protoPat[[full.var]] <- c(0,(rowSums(outer(rep(1, nrow(protoPat)),
          pattern == "full")*protoPat$LHartlDynf)/max_per_int)[-1])
    } else {
        CI   <- cumsum(pattern == "low")
        CII  <- cumsum(pattern == "mid")
        CIII <- cumsum(pattern == "full")
        protoPat[[low.var]]  <- lag(c(CI, rep(CI[maxdays.tdc],
          nrow(protoPat)-maxdays.tdc)), default = 0)
        protoPat[[mid.var]]  <- lag(c(CII, rep(CII[maxdays.tdc],
          nrow(protoPat)-maxdays.tdc)), default= 0)
        protoPat[[full.var]]  <- lag(c(CIII, rep(CIII[maxdays.tdc],
          nrow(protoPat)-maxdays.tdc)), default = 0)
      }
  }
  protoPat
}

#' Perform p-value calculation for difference between two hypothetical patients.
#'
#' @param effectname Name of variable(s) for which the effect should be extracted.
#' \code{"AdequacyCals"}.
#' @import mgcv
#' @importFrom stats predict
#' @keywords internal
#' @export
testEffectDiffs <- function(
  m,
  newdata1,
  newdata2,
  effectname = "AdequacyCals",
  debug      = FALSE,
  mod.df     = FALSE) {

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
  dropcols <- which(apply(X, 2, function(x) all(x == 0)))
  if(length(dropcols) > 0){
    useColumns <- useColumns[-dropcols]
    X <- X[,-dropcols]
    covCoefs <- covCoefs[-dropcols, -dropcols]
  }


  df <- sum(m$edf[useColumns])
  if (!is.null(m$edf1)) df <- sum(m$edf1[useColumns])
  df <- min(ncol(X), df)
  if(mod.df) df <- df/2
  if (m$scale.estimated) {
    rdf <-  length(m$y) - sum(m$edf)
  } else rdf <- -1

  if(!debug) {
    mgcv:::testStat(p = m$coefficients[useColumns], X = X, V = covCoefs,
      rank = df, type=0, res.df=rdf)
  } else {
    list(p = m$coefficients[useColumns], X = X, V = covCoefs,
      rank = df, type=0, res.df=rdf)
  }

}

#' The modus of a categorical variable
#'
#' @keywords internal
#' @export
modus <- function(f){
  freqs <- table(f)
  names(freqs)[which.max(freqs)]
}

#' Extract median covariate information from columns
#'
#' @importFrom stats median
#' @importFrom stats coefficients
#' @keywords internal
#' @export
median_patient <- function(
    data,
    patient.vars  = c("ApacheIIScore", "Age", "BMI", "Year", "DiagID2",
        "AdmCatID", "Gender", "freqInMVd2to4", "freqPFd2to4", "freqOId2to4",
        "freqPNd2to4"),
    random.effect = "CombinedicuID",
    re.byvar      = "icuByDummy") {

    median.x <- data.frame(
        lapply(
            subset(data, select=patient.vars),
            function(x) {
                if(is.factor(x)){
                    modus(x)
                } else {
                    median(x, na.rm=TRUE)
                }
            }))
    if(!is.null(random.effect)) {
    # pic any RE factor level and set according by variable to 0
        median.x[[random.effect]] <- factor(
            x = levels(data[[random.effect]])[1])
        median.x[[re.byvar]] <- 0
    }

    # return
    median.x

}
