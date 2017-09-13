## in order to simulate data we need a "full" data set, i.e. a data set
# in which all intervals are observed for each patient in the data

m <- readRDS("../../runModelBatchJobs/gamBatch/results/fullExpertt1Calories.Rds")
orig.data <- readRDS("../../dataCurrent/nutriOrigSmall.Rds")

vars.to.update = c("CombinedID", "intmid", "offset", "PatientDied",
	"event", "Survdays", "intMat")
brks <- c(0:30)
ints <- c(4:30)

comp.data <- complete_data(m, orig.data, vars.to.update=vars.to.update,
	brks=brks, ints=ints)


dir.create("input", showWarnings=FALSE)

saveRDS(comp.data, "input/completeData.Rds")