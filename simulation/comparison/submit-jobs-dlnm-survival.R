## load registry
library(batchtools)
reg <- loadRegistry("dlnm-tv-surv-registry/", writeable = TRUE, work.dir = getwd())

## set multicore or other parallel options
reg$cluster.functions  <- makeClusterFunctionsMulticore(ncpus=20)

#### DLNM for survival data
exp_dlnm_ped <- findExperiments(
	prob.name = "dlnm_sim_ped",
	algo.name = "pam_dlnm_ped")

# t1_pam_ped <- testJob(id=exp_dlnm_ped[1,1])
submitJobs(ids = findNotDone(exp_dlnm_ped[,1]))
waitForJobs()

#### DLNM time-varying PED
reg$cluster.functions  <- makeClusterFunctionsMulticore(ncpus=20)
exp_dlnm_pam_tv <- findExperiments(
	prob.name = "sim_dlnm_ped_tv",
	algo.name = "pam_dlnm_ped_tv")

# t1_pam_ped_tv <- testJob(id=exp_dlnm_pam_tv[1,1])
submitJobs(ids = findNotDone(exp_dlnm_pam_tv[,1]))
waitForJobs()