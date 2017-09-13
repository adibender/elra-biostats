##' A function that generates a simulated data set.
##' @static A list of object needed to simulate new data sets.
##' @param Name of the model from which hazards for complete data were obtained.
##' These hazards are used to simulate new survival times. Data set is then
##' adjusted according to simulated survival times.
##' @return \code{dynamic} part of simulation.

simPED_wrapper <- function(
  static, 
  haz.name, 
  brks=0:30,
  ints=4:30,
  vars.to.update=c("CombinedID", "intmid", "offset", "PatientDied",
    "event", "Survdays", "intMat"),
  ...) {

  dynamic <- sim_ped(
    static[["hazards.list"]][[haz.name]],
    static[["complete.data"]], 
    brks=brks, 
    ints=ints, 
    vars.to.update=vars.to.update, ...)

  # return

  dynamic[, c("CombinedID", "event", "intMat", "DaysMat", "intmid", "offset", 
    "AdequacyCalsTot30To70", "AdequacyCalsTotAbove70", "LHartlDynf", "LSimple460f",
    "ApacheIIScore")]

}


simPED_wrapper2 <- function(static, ...) {

  sim_ped(static[["hazards"]], static[["complete.data"]], ...)

}


## for evaluation of dependency of empirical rejection rate on df of smooth
gen_simple_data <- function(beta0=1, beta1=2, n=1000) {

  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y <- beta0 + beta1*x1 + rnorm(n)

  data.frame(y=y, x1=x1, x2=x2)

}

gen_simple_data_re <- function(beta0=1, beta1=2, n=1000, ngroup=50, sd.re=2) {

  id <- factor(rep(seq_len(ngroup), each=n/ngroup))
  ranef <- rep(rnorm(ngroup, 0, sd.re), each=n/ngroup)

  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y <- beta0 + beta1*x1 + ranef + rnorm(n)

  data.frame(y=y, x1=x1, x2=x2, id)

}


get_hazards_complete <- function(
  form = event ~ 
    te(intMat, DaysMat, by = I(AdequacyCalsTotAbove70 * LHartlDynf), 
      bs = "ps", m = list(c(2, 2), c(2, 2)), id = "cal") + 
    te(intMat, DaysMat, by = I(AdequacyCalsTot30To70 * LHartlDynf), 
      bs = "ps", m = list(c(2, 2), c(2, 2)), id = "cal"),
  ped, 
  ped.complete) {

  library(mgcv)

  m <- gam(form, family=poisson(), offset=offset, method="REML", data=ped, 
    control=gam.control(trace=TRUE, maxit=20))

  hazards <- predict(m, newdata=ped.complete, type="response")

  return(hazards)

}
