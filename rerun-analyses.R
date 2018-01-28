## title: Reproduce analyses in Bender et. al (2017)

# Note, if you just want to check if the code runs, without wanting to replicate
# the exact results, set n_simA and n_simB to lower numbers below.
# n_simA = 50 and n_simB = 100 should be enough for the graphs to look very
# similar to the ones in the publication


devtools::install("elrapack/")
# devtools::install_github("dgrtwo/gganimate")
devtools::install_github("adibender/pammtools@v0.0.3.2")
library(elrapack)
library(parallel)

## data generation

# (assuming your working directory is set to the main folder):
setwd("dataGenerationScripts/")
source("dataImportFromSAStoCleaned.R", echo=TRUE)
rm(list=ls())
gc()

## Model estimation

# WARNING This part can be especially memory intensive (30GB RAM or more per model)
# Probably best to run on a server
setwd("../runModelBatchJobs/")
source("combinations.R", echo=TRUE)
source("createSubmitModelFiles.R", echo=TRUE)
source("submitModelsGAM.R"       , echo=TRUE)
source("submitModelsGAMHosp.R"   , echo=TRUE)
rm(list=ls())
gc()

## Simulation Part B
# preparation
setwd("../")
source("protocols.R", echo=TRUE)
setwd("simulation/modelEvaluation/")
source("createCompleteData.R", echo=TRUE)
# The following code fits 3 models to obtain parameters for simulation
# each can take up to 25GB of memory
source("createSettingParameters.R", echo=TRUE)
# Actual simulation: # each job takes up to 5GB  of RAM
# To reduce the run-time, set a lower number of replications (per sub-setting)
n_simB <- 20 # (replication per setting, of which there are 16)
# set to value below for full replication
# n_simB <- 500
source("setupBatchExperimentHypocalorics.R", echo=TRUE)
source("submitJobs.R", echo=TRUE)
source("runSimEval.R", echo=TRUE)

rm(list=ls())
gc()

## Simulation Part A

# Scenario (1): each Job takes up to 5GB of RAM
# Scenario (2): each Job takes up to 10GB of RAM

setwd("../comparison/")
source("create-static-dlnm-survival.R", echo=TRUE)
# To reduce run-time, set a lower number of replications (per scenario)
n_simA <- 20
# set to value below for full replication
# n_simA <- 500
source("setup-batch-dlnm-survival.R", echo=TRUE)
source("submit-jobs-dlnm-survival.R", echo=TRUE)
rm(list=ls())
gc()

## Recreate Graphs and Tables for publication

# (assumes that previous steps ran without errors)
setwd("../../paper/lda")
source("prepLagLeadWindow.R", echo=TRUE)
source("prepModelResults.R" , echo=TRUE)
source("prepSimStats.R"     , echo=TRUE)
source("prepPvals.R"        , echo=TRUE)
source("simFunctionShape.R" , echo=TRUE)


## SessionInfo
sessionInfo()


## SessionInfo2
devtools::session_info()
