library(BatchJobs)


#### setup basic variables
# define project path
# project.path <- "/media/Data/subversion/hypocalorics/trunk/"
project.path <- "../"
# path to data
XXdatapathXX #
data.path <- paste0(project.path, data.folder, "/")

## setup formulae
# base formula (confounders)
XXbaseformulaXX


## nutrition formula
nutri.formula <- '+ te(intMat, DaysMat, by = I(AdequacyCalsTotAbove70 * LHartlDynf),
    bs = "ps", m=list(c(2,1),c(2,1)), id = "cal") +
  te(intMat, DaysMat, by = I(AdequacyCalsTot30To70 * LHartlDynf),
    bs = "ps", m=list(c(2,1),c(2,1)), id = "cal")'
nutri.formula2 <- '+ te(intMat, DaysMat, by = I(AdequacyCalsTotAbove70 * LHartlDynf),
    bs = "ps", m=list(c(2,2),c(2,2)), id = "cal") +
  te(intMat, DaysMat, by = I(AdequacyCalsTot30To70 * LHartlDynf),
    bs = "ps", m=list(c(2,2),c(2,2)), id = "cal")'

nutri.formulae <- list(
  t1    = nutri.formula,
  t2    = nutri.formula2)


data.files  <- list(
  full             = paste0(data.path, "nutriOrigSmall.Rds")
)

# load different subgroups/model specifications for which a model should be fit
combinations <- readRDS("combinations.Rds")
if(grepl(".*Hosp", data.folder)) {
  combinations <- combinations["fullExpertt1Calories", , drop=FALSE]
}



## setup project specifics
# load server configuration (for batch jobs )
loadConfig(paste0(project.path, "BatchJobsLocal.R"))
## set paths and environment variables
XXregnameXX # reg.name <- "bam182"
XXrunmodel.functionXX #
function.path <- paste0(project.path, "runModelBatchJobs/", function.name)

## unlink in case study is rerun
unlink(reg.name, recursive = TRUE)
reg <- makeRegistry(
  id        = reg.name,
  file.dir  = reg.name,
  src.files = function.path)

save.path <- paste0(reg$file.dir, "/results/")
dir.create(save.path, showWarnings = FALSE)

fun.args <- list(
  all.combinations = combinations,
  base.formula     = base.formula,
  nutri.formulae   = nutri.formulae,
  data.files       = data.files,
  save.path        = save.path)

## runModel in runOneModelGAM.R
batchMap(reg, runModel, rownames(combinations), more.args = fun.args,
  use.names = TRUE)

# submit jobs and wait for them to finish
submitJobs(reg)
waitForJobs(reg)
