library(checkmate)
library(mgcv)

## path to models on full data set/sensitivity analyses/subgroup analyses
res.path <- "../../runModelBatchJobs/gamBatch/results"
assert_directory(res.path, "r")

## file names for the various models (needed later to load the models and
# extract the coefficients)
mod.files <- list.files(
    path       = res.path,
    full.names = TRUE)
# complete data needed to simulate hazards for alle subjects and all intervals
complete.data <- readRDS("input/completeData.Rds")
# nutri data used to obtain coefficients for reduced models
nutri <- readRDS("../../dataCurrent/nutriOrigSmall.Rds")

## Define model specifications that will be used to simulate data and
## that will be fit to simulated data
# (currently only Expertt1, Expertt2, Statict1 and NULL relevant)
form.Expertt1 <- event ~ s(intmid, bs="ps") + te(ApacheIIScore, intmid) +
    te(intMat, DaysMat, by=I(AdequacyCalsTot30To70*LHartlDynf),
        bs="ps", m=list(c(2, 1), c(2, 1)), id="cal") +
    te(intMat, DaysMat, by=I(AdequacyCalsTotAbove70*LHartlDynf),
        bs="ps", m=list(c(2, 1), c(2, 1)), id="cal")
form.Statict1 <- event ~ s(intmid, bs="ps") + te(ApacheIIScore, intmid) +
    te(intMat, DaysMat, by=I(AdequacyCalsTot30To70*LSimple460f),
        bs="ps", m=list(c(2, 1), c(2, 1)), id="cal") +
    te(intMat, DaysMat, by=I(AdequacyCalsTotAbove70*LSimple460f),
        bs="ps", m=list(c(2, 1), c(2, 1)), id="cal")
form.Expertt2 <- event ~ s(intmid, bs="ps") + te(ApacheIIScore, intmid) +
    te(intMat, DaysMat, by=I(AdequacyCalsTot30To70*LHartlDynf),
        bs="ps", m=list(c(2, 2), c(2, 2)), id="cal") +
    te(intMat, DaysMat, by=I(AdequacyCalsTotAbove70*LHartlDynf),
        bs="ps", m=list(c(2, 2), c(2, 2)), id="cal")
form.Statict2 <- event ~ s(intmid, bs="ps") + te(ApacheIIScore, intmid) +
    te(intMat, DaysMat, by=I(AdequacyCalsTot30To70*LSimple460f),
        bs="ps", m=list(c(2, 2), c(2, 2)), id="cal") +
    te(intMat, DaysMat, by=I(AdequacyCalsTotAbove70*LSimple460f),
        bs="ps", m=list(c(2, 2), c(2, 2)), id="cal")
form.null <- event ~ s(intmid, bs="ps") + te(ApacheIIScore, intmid)

form.list <- list(
    Expertt1     = form.Expertt1,
    Statict1     = form.Statict1,
    Expertt2     = form.Expertt2,
    Statict2     = form.Statict2,
    NoneNoneNone = form.null)

# for each model (one list entry), save vector of hazards, each hazard for
# patient i, interval j
# need to be evaluated iteratively, at leas parallel:::mclapply throughs errors

## get coeffiicients for adequacy vars for simulation from models on real data
m.fullExpertt1Calories <- readRDS(mod.files[1])
m.fullExpertt2Calories <- readRDS(mod.files[2])
m.fullL460t1Calories   <- readRDS(mod.files[3])

ind.adeq.full <- grep("Adequacy", names(coef(m.fullExpertt1Calories)))
coef.dyn1  <- coef(m.fullExpertt1Calories)[ind.adeq.full]
coef.dyn1  <- c(coef.dyn1[26:50], 1.2*coef.dyn1[26:50])
coef.dyn2  <- coef(m.fullExpertt2Calories)[ind.adeq.full]
coef.dyn2  <- c(coef.dyn2[26:50], 1.2*coef.dyn2[26:50])
coef.stat1 <- coef(m.fullL460t1Calories)[ind.adeq.full]
coef.stat1 <- c(coef.stat1[26:50], 1.2*coef.stat1[26:50])

rm(m.fullExpertt1Calories, m.fullExpertt2Calories, m.fullL460t1Calories)
gc()

# subset, b/c for time-constant covars first interval contains all information
med.patient <- median_patient(subset(complete.data, intmid==4.5))

## create compare list
maxdays.tdc <- ncol(complete.data[["DaysMat"]])
# can pick any ID, since in complete data, all have same number of observations
maxn.intervals <- length(unique(complete.data[["intmid"]]))
# define protocols which should be used in comparisons
protocols.df <- readRDS("../../protocolsDF.Rds")

protocols.to.compare <- readRDS("../../protocolsToCompare.Rds")


hazards.list <- mclapply(seq_along(form.list), function(z) {

  print(form.list[[z]])

  m <- gam(
    formula = form.list[[z]],
    data    = nutri,
    offset  = offset,
    family  = poisson(),
    method  = "REML")
  cm <- coef(m)
  cm[1] <- -4.5
  ind.adeq <- grep("Adequacy", names(cm))
  if(names(form.list)[z] == "Expertt1")
    cm[ind.adeq] <- coef.dyn1
  if(names(form.list)[z] == "Expertt2")
    cm[ind.adeq] <- coef.dyn2
  if(names(form.list)[z] == "Statict1")
    cm[ind.adeq] <- coef.stat1

  X <- predict(m, newdata=complete.data, type="lpmatrix")
  Vp <- m$Vp


  sub.list  <- list(
    hazards = as.numeric(exp(X %*% cm)),
    formula = form.list[[z]],
    coefs   = cm,
    Vp      = Vp,
    lp.list = lapply(protocols.df, function(pat) {
      proto_lp(
        pattern  = pat,
        model    = m,
        ped      = complete.data,
        median.x = med.patient)
    }))

  return(sub.list)

}, mc.cores=5)


formulas.list <- lapply(hazards.list, "[[", i="formula")
# need to do this, otherwise saving and reading takes forever
formulas.list <- lapply(formulas.list, function(z) {
  attr(z, ".Environment") <- NULL
  z
})
coefs.list    <- lapply(hazards.list, "[[", "coefs")
Vp.list       <- lapply(hazards.list, "[[", "Vp")
Xcompare.list <- lapply(hazards.list, "[[", "lp.list")
hazards.list  <- lapply(hazards.list, "[[", "hazards")

names(formulas.list) <- names(coefs.list) <- names(hazards.list) <-
    names(Vp.list) <- names(Xcompare.list) <- names(form.list)

saveRDS(Xcompare.list , "input/XcompareList.Rds")
saveRDS(formulas.list , "input/formulasList.Rds")
saveRDS(coefs.list    , "input/coefsList.Rds")
saveRDS(Vp.list       , "input/VpList.Rds")
saveRDS(hazards.list  , "input/hazardsList.Rds")

## final object, provided to addProblem.
static <- list(
    complete.data        = complete.data,
    hazards.list         = hazards.list,
    formulas.list        = formulas.list,
    coefs.list           = coefs.list,
    Vp.list              = Vp.list,
    Xcompare.list        = Xcompare.list,
    protocols.to.compare = protocols.to.compare,
    protocols.df         = protocols.df,
    maxdays.tdc          = ncol(complete.data[["DaysMat"]]),
    intmid               = sort(unique(complete.data$intmid)))


saveRDS(static, "input/static.Rds")