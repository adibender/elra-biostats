# @Author: andreas.bender@stat.uni-muenchen.de
# @Date:   2015-02-04 18:14:03
# @Last Modified by:   andreas.bender@stat.uni-muenchen.de
# @Last Modified time: 2015-06-10 15:39:49




#' Select the interval in which given t falls and return according data set, id
#' or index

get_tid <- function(
  df,
  t,
  return.val = "df",
  id.var     = "CombinedID",
  start.var  = "tstart",
  end.var    = "tend",
  select     = NULL) {

  diff.tend <- sign(df[[end.var]] - t)
  diff.start <- sign(t - df[[start.var]])
  ind <- which((diff.start == 1) & (diff.tend != -1))

  if (return.val == "df") {
    if (is.null(select)) df[ind, ]
    else df[ind, select]
  } else if(return.val == "id") {
    df[[id.var]][ind]
  } else{
    ind
  }

}

#' Given data set in long format and ordered by an id, returns last
# observationsfor according id

last_ind <- function(df, id.var = "CombinedID") {

  cumsum(rle(df[[id.var]])$lengths)

}

#' Given data set in long format, returns data set of last observations for an
# id and only for timepoints before a specified time t.
get_iddf <- function(df, t, id.var = "CombinedID") {

  df <- subset(df, tstart < t)
  df[last_ind(df, id.var = id.var), ]

}

#' Returns data frame containing predicted hazards and selected variables
# up until a specified time t.
get_hazard <- function(
  object,
  newdata,
  t,
  start.var = "tstart",
  vars.keep = c("CombinedID", "tstart", "tend", "survtime")) {

  risk.set <- subset(newdata, newdata[[start.var]] < t)

  # get predictions (hazard rates) for all time points before 't' of interest
  mu.ij <- predict(object, newdata = risk.set, type = "response")

  risk.set$lambda.ij <- mu.ij
  risk.set$Lambda.ij <- unlist(with(risk.set,
    sapply(split(lambda.ij, CombinedID), cumsum)))
  risk.set$St.ij <- exp(-risk.set$Lambda.ij)

  hazards.mat <- risk.set[, c(vars.keep, "lambda.ij", "Lambda.ij", "St.ij")]
  class(hazards.mat) <- c("hazards.mat", class(hazards.mat))

  hazards.mat

}

#' Transforms a data frame containing Survival probabilities in long format to
# wide format
get_survprob <- function(St.df) {

  library(reshape2)
  dcast(St.df, CombinedID~ tend, value.var = "St.ij")

}

#' Tests if provided survival probabilities are monotnously decreasing.
test_monotony <- function(survprobs) {

  survprobs <- na.omit(survprobs)
  all(survprobs == cummin(survprobs))

}

#' Given a matrix of survival probabilities (rows = ids, cols = times), tests
# if all rows contain monotonously decreasing survival probabilities.

test_survprob <- function(pmat, t) {

  if(all(is.na(pmat))) stop("All entries NA")
  rows <- apply(pmat, 1, test_monotony)
  all(rows)

}

##NOTE: as of now pmat as returned by get_hazard/get_survprob contains NA for
# intervals that person didn't live long enough to see... altough
# estimation/prediction of S(t) for these intervals would be theoretically
# possible (as in simulation studies). This would invoke extending the
# data sets of individudal patients (partially imputing nutrition protocols
# would be required)
# Since these imputations are no performed as of now, data provied to cindex
# function can only contain patients who survived until the begining of the
# interval for which cindex should be calculated (selection by get_tid)

#' Calculate C-Index at time t
#'
#' For specified time t and given survival probabilities for this time point,
#' retuns (apparent) C-index evaluated at the respective time t.
#'
#' @param pmat Matrix of surival probabilities (rows = patients, cols = time
# points).
#' @param newdata For the apparent C-Index same data that was used to calculate
#' \code{pmat}.
#' @param t Time point at which C-Index should be evaluated.
#' @param form.ipcw Formula for the model which will be used to calculate inverse
#' probabilities of censoring weights.
#' @param cens model If "none", no IPCW correction is performed.
cindex_t <- function(
  pmat,
  newdata,
  t,
  form.ipcw  = Hist(survtime, event) ~ 1,
  cens.model = "none",
  ...) {

  library(pec)
  library(survival)
  library(prodlim)

  if(!is.matrix(pmat)) pmat <- as.matrix(pmat)
  pmat <- na.omit(pmat)

  cindex(
    object     = pmat,
    data       = get_tid(newdata, t),
    eval.times = t,
    # pred.times = t,
    formula    = form.ipcw,
    cens.model = cens.model)$AppCindex[[1]]

}

#' Returns (apparent) C-Index for each time point specified in \code{times}.

#' @param object A matrix with survival probabilities (rows = patients, cols =
#' time).
#' @param newdata Data set on which C-Index will be evaluated.
#' @param times Times at which C-Index should be evaluated.
#' @param form.ipcw Formula for the model which will be used to calculate inverse
#' probabilities of censoring weights.
#' @param cens model If "none", no IPCW correction is performed.
#' @rdname cindex
cindex_times <- function(
  object,
  newdata,
  times,
  form.ipcw  = Hist(survtime, event) ~ 1,
  cens.model = "none", ...) {

  cind <- vapply(
    times,
    function(t) {
      cindex_t(
        pmat       = object[, colnames(object) == t],
        newdata    =  newdata,
        t          = t,
        form.ipcw  = form.ipcw,
        cens.model = cens.model, ...)},
    FUN.VALUE=numeric(1))

  data.frame(cindex = cind, times = times)

}

#' @inheritParams cindex_times
#' @rdname cindex
#' @export
cindex_wrapper <- function(
  model,
  data,
  times = NULL,
  form.ipcw = Hist(survtime, event) ~ 1,
  cens.model = "none",
  tvar = "tend", ...) {


  if(is.null(times)) {
    if(!(tvar %in% colnames(data))) stop(paste0(tvar, " not in data set!"))
    times <- unique(data[[tvar]])
  }

  pmat <- get_survprob(get_hazard(model, data, max(times)))[, -1]
  cindex_times(pmat, data, times, form.ipcw = form.ipcw, cens.model = cens.model, ...)

}
