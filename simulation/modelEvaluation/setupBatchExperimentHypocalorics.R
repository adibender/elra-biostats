library(BatchExperiments)

# setup cluster/server configuration
#loadConfig("../../BatchJobs.R")
loadConfig("../../BatchJobsLocal.R")

# load inputs
static <- readRDS("input/static.Rds")

## create registry for evaluation of gam models for ELRA effects
# source.dirs contains functions that will aid simulation, e.g. simPED takes
# is the "dynamic" part in BatchExperiments jargon
reg <- makeExperimentRegistry(
    id        = "modelEvalELRA",
    src.files = c("algorithms.R", "problems.R"),
    packages  = c("mgcv", "elrapack"))


########### Problem 1:
## fit all models to data generated from fullExpertt1Calories model

## add algoirthm (bam_wrapper defined in algorithms.R script)
source("algorithms.R")
addAlgorithm(reg, id="bam", fun=bam_wrapper, overwrite=TRUE)
addAlgorithm(reg, id="gam", fun=gam_wrapper, overwrite=TRUE)

## add problem (dynamic function defined in problems.R script)
source("problems.R")
addProblem(
    reg,
    id        = "dataFromExpertt1",
    static    = static,
    dynamic   = simPED_wrapper,
    seed      = 29082016,
    overwrite = TRUE)

# create designs for Problem 1
dataFromExpertt1 <- makeDesign(
    id         = "dataFromExpertt1",
    exhaustive = list(haz.name=c("Expertt1")))

gam.design.dataFromExpertt1 <- makeDesign(
    id         = "gam",
    exhaustive = list(form=names(static$formulas.list)[c(1, 2, 3, 4)]))

# add Experiment for Problem 1
addExperiments(
    reg,
   repls         = n_simA,
    prob.design  = dataFromExpertt1,
    algo.design  = list(gam.design.dataFromExpertt1),
    skip.defined = TRUE)


########### Problem 2:
## - Data is generated from model with static lag lead window

## add Problem
addProblem(
    reg,
    id        = "dataFromStatict1",
    static    = static,
    dynamic   = simPED_wrapper,
    seed      = 20150910,
    overwrite = TRUE)

# make design for Problem 2
dataFromL460t1.design <- makeDesign(
    id = "dataFromStatict1",
    exhaustive = list(haz.name="Statict1"))

gam.design.dataFromL460t1 <- makeDesign(
    id         = "gam",
    exhaustive = list(form = names(static$formulas.list)[c(1, 2, 3, 4)]))

# ad Experiments for Problem 2:
addExperiments(
    reg,
   repls         = n_simA,
    prob.design  = dataFromL460t1.design,
    algo.design  = gam.design.dataFromL460t1,
    skip.defined = TRUE)


########### Problem 3:
## - Data is generated from fullExpertt2Calories

## add Problem
addProblem(
    reg,
    id        = "dataFromExpertt2",
    static    = static,
    dynamic   = simPED_wrapper,
    seed      = 20150912,
    overwrite = TRUE)

# make design for Problem 3
dataFromExpertt2.design <- makeDesign(
    id = "dataFromExpertt2",
    exhaustive = list(haz.name=c( "Expertt2")))

gam.design.dataFromExpertt2 <- makeDesign(
    id         = "gam",
    exhaustive = list(form = names(static$formulas.list)[c(1, 2, 3, 4)]))
addExperiments(
    reg,
   repls         = n_simA,
    prob.design  = dataFromExpertt2.design,
    algo.design  = gam.design.dataFromExpertt2,
    skip.defined = TRUE)


########### Problem 4:
## - Data is generated from fullNoneNoneNone model,
# i.e. model without terms associated with nutrition variables
# should simulate setting where coefs vor nutrition effects = 0

## add Problem
addProblem(
    reg,
    id        = "dataFromNoneNoneNone",
    static    = static,
    dynamic   = simPED_wrapper,
    seed      = 20150914,
    overwrite = TRUE)

# make design for Problem 3
dataFromNoneNoneNone.design <- makeDesign(
    id = "dataFromNoneNoneNone",
    exhaustive = list(haz.name="NoneNoneNone"))

gam.design.dataFromNoneNoneNone <- makeDesign(
    id         = "gam",
    exhaustive = list(
        form = c(
            "Expertt1",
            "Statict1",
            "Expertt2",
            "Statict2")))
addExperiments(
    reg,
   repls         = n_simA,
    prob.design  = dataFromNoneNoneNone.design,
    algo.design  = gam.design.dataFromNoneNoneNone,
    skip.defined = TRUE)
