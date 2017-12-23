## load registry
library(batchtools)
reg <- loadRegistry("dlnm-tv-surv-registry/", writeable = TRUE, work.dir = getwd())

## set multicore or other parallel options
reg$cluster.functions  <- makeClusterFunctionsMulticore(ncpus = 20)

#### PAM WCE for DLNM survival data
exp_wce_ped <- findExperiments(
  prob.name = "dlnm_sim_ped",
  algo.name = "pam_wce_ped")

t1_pam_wce_dlnm <- testJob(id = exp_wce_ped[1,1])
submitJobs(ids = findNotDone(exp_wce_ped[,1]))
waitForJobs()

#### PAM DLNM for DLNM survival data
exp_dlnm_ped <- findExperiments(
	prob.name = "dlnm_sim_ped",
	algo.name = "pam_dlnm_ped")

t1_pam_dlnm_dlnm <- testJob(id=exp_dlnm_ped[1,1])
submitJobs(ids = findNotDone(exp_dlnm_ped[, 1]))
waitForJobs()


#### PAM DLNM fit on TV DLNM
reg$cluster.functions  <- makeClusterFunctionsMulticore(ncpus=30)
exp_pam_dlnm_tvdlnm <- findExperiments(
  prob.name = "sim_dlnm_ped_tv",
  algo.name = "pam_dlnm_ped")

t1_pam_dlnm_tvdlnm <- testJob(id=exp_pam_dlnm_tvdlnm[1,1])
submitJobs(ids = findNotDone(exp_pam_dlnm_tvdlnm[, 1]))
waitForJobs()

#### PAM TV DLNM fit on TV DLNM
reg$cluster.functions  <- makeClusterFunctionsMulticore(ncpus=30)
exp_dlnm_pam_tv <- findExperiments(
	prob.name = "sim_dlnm_ped_tv",
	algo.name = "pam_dlnm_ped_tv")

t1_pam_ped_tv <- testJob(id=exp_dlnm_pam_tv[1,1])
submitJobs(ids = findNotDone(exp_dlnm_pam_tv[, 1]))
waitForJobs()
