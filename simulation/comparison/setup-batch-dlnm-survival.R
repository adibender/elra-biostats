library(batchtools)

reg <- makeExperimentRegistry(
	"dlnm-tv-surv-registry",
	packages = c("mgcv", "magrittr", "dplyr", "purrr", "pam", "dlnm", "ggplot2"),
	source   = c("problems-dlnm-survival.R", "algorithms-dlnm-survival.R"),
	seed     = 10072017)

reg <- loadRegistry("dlnm-tv-surv-registry/", writeable=TRUE)

reg$cluster.functions = makeClusterFunctionsMulticore(ncpus=20)
saveRegistry()

#### DLNM for survival data
addProblem(
	name = "dlnm_sim_ped",
	data = readRDS("input/static_dlnm_ped.Rds"),
	fun  = sim_dlnm_ped)

addAlgorithm(
	name = "pam_dlnm_ped",
	fun = pam_dlnm_ped)

addExperiments(
	prob.designs = list(
		dlnm_sim_ped = data.frame()),
	algo.designs = list(
		pam_dlnm_ped = data.frame(
			debug = FALSE)),
	repls = n_simB)


#### time-varying DLNM for survival data
addProblem(
	name = "sim_dlnm_ped_tv",
	data = readRDS("input/static_dlnm_ped_tv.Rds"),
	fun = sim_dlnm_ped_tv)
addAlgorithm(
	name = "pam_dlnm_ped_tv",
	fun = pam_dlnm_ped_tv)

addExperiments(
	prob.designs = list(
		sim_dlnm_ped_tv = data.frame()),
	algo.designs = list(
		pam_dlnm_ped_tv = data.frame(
			debug = FALSE)),
	repls = n_simB)