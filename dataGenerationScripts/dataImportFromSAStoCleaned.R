## Original data was obatained as a SAS data base together with a SAS script
# that exported the data to 3 files in csv format:
# - ICU.csv: Contains information on individual ICUs, e.g. country and other ICU
# specific variables
# - daily.csv: Contains daily information regarding the feeding protocol (at
# most 12 # observations (study days 1-12) per patient, as well mechanical
# ventilation, etc. )
# - patient.csv: Contains patient related data of not time-dependent nature
# (i.e. patient characteristics known at time of admission like weight, age,
# etc.)

## The following steps are necessary to create data sets that were used for
# final data analysis
#
# 1) importHypocaloric.R: imports data from "*.csv" to R and performs
# initial data cleaning etc. The processed data sets are stored
# in "*.Rds" files ICU, patient, daily and mergedAndCleanedData
# which is merely a merge of the other 3 data sets (not included due to
# anonymization. The data sets included in the "data" folder resulted from
# this original import script)
#
# 2) createPEMdata.R: transforms the data into PED (Piece-wise Exponential Data)
# format, i.e. one row for each observed interval per patient. Specification of
# the intervals (interval lengths) can be passed to the according
# transformation functions. Time-dependent covariates that are used later are
# also specified here.
# Output is stored in "dataPEM.Rds"
#
# 3) selectValidCases.R: performs validity checks/applies exclustion criteria
# on the newly created PED and returns a subset with valid cases "validCases.Rds"
#
# 4) adjustRawDataToValidCases.R: removes invalid cases as computed in the
# previous step, from the raw data sets (patient, daily, ICU,  mergedAndCleanedData),
# creates data sets patientVC, dailyVC, etc.
#
# 5) createLagLeadNutritionVars.R: adds nutrition variables to
# "validCases" data. Different definitions are used to enable subsequent sensitivty
# analyses.
# Output is stored in "dataWithNutrition.Rds".
#
# 6) prepData.R performs some additional data preparation
#
# 7) createSubsetsForModels.R: subsets the full data set
# (which speeds up data read-in process and reduces memory load).
#
# 7) checkDimensions.R: A final run through the created data sets to check
# whether dimensions (rows and columns) correspond to the ones obtained by
# the authors. If not, an error will be thrown pointing to the data sets with
# conflicting dimensions. If this happens, it means that something went wrong
# (maybe differing package versions).


## load packages with additional functionality
library(lubridate) # handling date formats/variables
library(checkmate) # checking inputs
library(dplyr)     # data manipulation library

# set project path
project.path <- "../"
# set path to files that should be sourced
path.to.functions     <- paste0(project.path, "functions/")
path.data.gen.scripts <- paste0(project.path, "dataGenerationScripts/")

## load helper functions
save.folder         <- "dataCurrent"
save.folder.hosp    <- "dataCurrentHosp"
save.path.data      <- paste0(project.path, save.folder, "/")
save.path.data.hosp <- paste0(project.path, save.folder.hosp, "/")
dir.create(save.path.data,      showWarnings = FALSE)
dir.create(save.path.data.hosp, showWarnings = FALSE)

# load basic setup
# defines some "global" variables like interval break points, etc.
source(paste0(project.path, "setup.R"))

## reading in data, preprocessing, sanity checks, new variables, merging data sets
# source(paste0(path.data.gen.scripts, "importHypocaloric.R"), echo=TRUE)

## create data in format suitable for piece-wise exponential models
source(paste0(path.data.gen.scripts, "createPEMdata.R"), echo=TRUE)

## apply exclusion criteria
source(paste0(path.data.gen.scripts, "selectValidCases.R"), echo=TRUE)

## remove excluded observations also from other data sets
source(paste0(path.data.gen.scripts, "adjustRawDataToValidCases.R"), echo=TRUE)

## create nutrition variables (need to be matrices) as well as lag/lead matrices
# that define which days of nutrition protocol can affect
source(paste0(path.data.gen.scripts, "createLagLeadNutritionVars.R"), echo=TRUE)

# some more adjustments, new variables needed for description
source(paste0(path.data.gen.scripts, "prepData.R"), echo=TRUE)

# create subset for sensitivity/subgroup analyses
# use only covariates that are needed for models to reduce RAM load when fitting
# models
source(paste0(path.data.gen.scripts, "createSubsetsForModels.R"), echo=TRUE)

## check if created data has expected dimensions
source(paste0(path.data.gen.scripts, "checkDimensions.R"), echo=TRUE)
