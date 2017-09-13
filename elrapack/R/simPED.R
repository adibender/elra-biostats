#' Simulate Data for simulation Part A
#'
#' Draws survival times from Piece-wise Exponential distribution and
#' return new data set in Piece-wise Exponential data (PED) format.
#'
#' @param hazards.ij Hazard for each subject and interval combination.
#' @param complete.data Full data set (i.e. with one row per interval)
#' for each subject and complete exposure history in PED format.
#' @import mgcv
#' @importFrom msm rpexp
#' @importFrom plyr colwise
#' @keywords internal
#' @export
sim_ped <- function(

    hazards.ij,
    complete.data,

    time.shift = 4,
    brks = c(0:12, seq(15, 55, by = 5), 61.1),
    ints = c(4:12, seq(15, 55,  by = 5), 61),

    vars.to.update = c("CombinedID", "intmid", "offset", "PatientDied",
        "event", "Survdays", "intMat"),
    vars.to.delete = c("inICU", "inMechVent"),
    re = TRUE,
    survvar="Survdays") {

    max.time <- max(ints)

    ## vars to keep
    adequacy.vars <- grep(".*Adequacy.*", colnames(complete.data), value = TRUE)
    L.vars <- grep("^L.*", colnames(complete.data), value = TRUE)
    vars.to.delete <- c(vars.to.update, vars.to.delete, adequacy.vars, L.vars)
    vars.to.keep <- colnames(complete.data)[!(colnames(complete.data) %in%
        vars.to.delete)]


    ## times
    time.int <- ints - time.shift

    ## get hazard rates and times
    valid.id <- complete.data$CombinedID
    ind.first <- which(!duplicated(valid.id))
    rle.id <- rle(valid.id)
    n.obs <- rle.id$lengths
    seq.id <- lapply(n.obs, seq_len)
    timez <- rep(time.int[-length(time.int)], times = length(n.obs))

    rt.l <- split(
        cbind.data.frame(
            rate = hazards.ij,
            t = timez),
        f = complete.data$CombinedID)
    rt.l <- lapply(rt.l, as.list)



    ## simulate new survival
    new.times <- time.shift + sapply(rt.l, do.call, what = rpexp)
    cens <- numeric(length(new.times))

    new.times[new.times > max.time] <- max.time
    cens <- new.times == max.time


    ## create new df and adjust offset + event
    n.obs.new <- findInterval(new.times, brks[-1:-time.shift], rightmost.closed=TRUE)
    seq.id.new <- lapply(n.obs.new, seq_len)
    id.new <- rep(unique(valid.id), times = n.obs.new)
    new.df <- data.frame(CombinedID = id.new, event = 0)

    ## fix adequacy to reflect new survival (locf if new.surv > old.surv)
    adequacy.new <- colwise(
        getAdequacyMat,
        ind.first = ind.first,
        n.obs.old = n.obs,
        n.obs.new = n.obs.new)(complete.data[, adequacy.vars])

    ## fix L vars
    L.new <- colwise(getLmat,
            n.obs.new = n.obs.new,
            n.obs.old = n.obs,
            ind.first = ind.first)(complete.data[, L.vars])

    ## fix DaysMat
    new.df$intMat <- getLmat(complete.data[, "intMat"], ind.first, n.obs, n.obs.new)


    ## first and last ind per patient (new data)
    ind.first.new <- which(!duplicated(new.df$CombinedID))
    ind.last.new <- cumsum(n.obs.new)

    ## insert intmid
    new.df$intmid <- getIntmid(brks[-1:-time.shift], seq.id.new)
    new.df$intmidF <- as.factor(new.df$intmid)

    ## insert offset
    new.df$offset <- getOffset(brks[-1:-time.shift], n.obs.new, seq.id.new,
            new.times)

    ## insert event
    new.df$event[ind.last.new] <- !cens * 1

    ## update patientDied and Survdays variables
    new.df$PatientDied <- rep(new.df$event[ind.last.new], times = n.obs.new)
    new.df[[survvar]] <- rep(new.times, times = n.obs.new)
    rep.ind <- rep(seq_along(n.obs.new), times = n.obs.new)
    new.df.constant <- complete.data[!duplicated(complete.data$CombinedID), vars.to.keep][rep.ind, ]

    if(!re) new.df.constant$icuByDummy <- 0
    ## return constant and updated parts of the data frame
    cbind(new.df, adequacy.new, L.new, new.df.constant)

}


#' @inherit sim_ped
getIntmid <- function(brks, seq.obs) {

    intmid <- (brks[-length(brks)] + brks[-1]) / 2
    unlist(lapply(seq.obs, function(z) intmid[z]))

}

#' @importFrom zoo na.locf
#' @inheritParams sim_ped
#' @rdname sim_ped
getAdequacyMat <- function(mat, ind.first, n.obs.old, n.obs.new) {

    fix.val   <- ifelse(is.numeric(mat), 0, FALSE)
    n.obs.old <- n.obs.old + 4
    n.obs.new <- n.obs.new + 4
    new.mat   <- mat[ind.first, ]
    seq.max   <- seq_len(ncol(mat))
    seq.new   <- lapply(n.obs.new, seq_len)
    seq.old   <- lapply(n.obs.old, seq_len)
    na.mat <- matrix(unlist(lapply(seq.new, "[", i = seq.max)),
            ncol = length(seq.max), byrow = TRUE)
    na.mat.old <- matrix(unlist(lapply(seq.old, "[", i = seq.max)),
            ncol = length(seq.max), byrow = TRUE)
    new.mat[is.na(na.mat.old)] <- NA
    new.mat <- t(na.locf(t(new.mat)))
    new.mat[is.na(na.mat)] <- fix.val

    new.mat[rep(seq_along(n.obs.new), times = n.obs.new - 4), ]

}

#' @inheritParams sim_ped
#' @rdname sim_ped
#' @keywords internal
#' @export
getLmat <- function(mat, ind.first, n.obs.old, n.obs.new) {

  ## here n.obs used in the sense of intervals survived in time 4:60
  max.int <- max(n.obs.old)
  L.max <- mat[ind.first[max.int]:(ind.first[max.int] + max.int - 1), ]
  seq.new <- unlist(lapply(n.obs.new, seq_len))

  L.max[seq.new, ]

}

#' @inheritParams sim_ped
#' @rdname sim_ped
#' @keywords internal
#' @export
getOffset <- function(brks, n.obs, seq.obs, surv) {

  ## create offset vec with according lengths
  off.set <- diff(brks)
  off.set <- lapply(seq.obs, function(z) off.set[z])

    ## compute offset in last interval
  offset.vec <- unlist(off.set)
  brks.vec <- unlist(lapply(seq.obs, function(z) brks[z]))
  last.ind <- cumsum(n.obs)
  offset.vec[last.ind] <- surv - brks.vec[last.ind]

  ## return offset
  log(offset.vec)

}

#' @inherit sim_ped
#' @rdname sim_ped
#' @keywords internal
#' @export
complete_data <- function(
  m, train,

  time.shift = 4,
  brks = c(0:12, seq(15, 55, by = 5), 61.1),
  ints = c(4:12, seq(15, 55,  by = 5), 61),

  vars.to.update = c("CombinedID", "intmid", "offset", "PatientDied",
    "event", "Survdays", "intMat"),
  vars.to.delete = c("inICU", "inMechVent"), hazards = NULL) {

  ## vars to keep
  adequacy.vars <- grep(".*Adequacy.*", colnames(train), value = TRUE)
  L.vars <- grep("^L.*", colnames(train), value = TRUE)
  vars.to.delete <- c(vars.to.update, vars.to.delete, adequacy.vars, L.vars)
  vars.to.keep <- colnames(train)[!(colnames(train) %in% vars.to.delete)]


    ## times
  time.int <- ints - time.shift
    #minus.ind <- -1:-time.shift (maybe, but reduces readability)
  if( !is.null(m$na.action) ) train <- train[-m$na.action, ]

    ## get hazard rates and times
  valid.id  <- train$CombinedID
  ind.first <- which(!duplicated(valid.id))
  rle.id    <- rle(valid.id)
  n.obs     <- rle.id$lengths
  seq.id    <- lapply(n.obs, seq_len)
  # if( is.null(hazards) ) hazards <- fitted(m)/exp(train$offset)


  ## create new df and adjust offset + event
  n.obs.new  <- rep(max(n.obs), times = length(n.obs))
  seq.id.new <- lapply(n.obs.new, seq_len)
  id.new     <- rep(unique(valid.id), times = n.obs.new)
  new.df     <- data.frame(CombinedID = id.new, event = 0)

  ## fix adequacy to reflect new survival (locf if new.surv > old.surv)
  adequacy.new <- colwise(getAdequacyMat, ind.first = ind.first,
    n.obs.old = n.obs, n.obs.new = n.obs.new)(train[, adequacy.vars])

  ## fix L vars
  L.new <- colwise(getLmat, n.obs.new = n.obs.new,
    n.obs.old = n.obs, ind.first = ind.first)(train[, grep("*.f", L.vars, value = TRUE)])
  # colnames(L.new) <- sub("f", "", colnames(L.new))

  ## fix DaysMat
  new.df$intMat <- getLmat(train[, "intMat"], ind.first, n.obs, n.obs.new)


  ## first and last ind per patient (new data)
  ind.first.new <- which(!duplicated(new.df$CombinedID))
  ind.last.new  <- cumsum(n.obs.new)

  ## insert intmid
  new.df$intmid  <- getIntmid(brks[-1:-time.shift], seq.id.new)
  new.df$intmidF <- as.factor(new.df$intmid)

  ## insert offset
  new.df$offset <- rep(ints[-1], times = length(n.obs.new)) -
  rep(ints[-length(ints)], times = length(n.obs.new))

  rep.ind <- rep(seq_along(n.obs.new), times = n.obs.new)
  new.df.constant <- train[!duplicated(train$CombinedID), vars.to.keep][rep.ind, ]

  ## return constant and updated parts of the data frame
  cbind(new.df, adequacy.new, L.new, new.df.constant)

}
