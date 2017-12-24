library(BatchExperiments)
# load configuration setup of parallel environment (server/cluster/etc.)

loadConfig("../../BatchJobsLocal.R")
#loadConfig("../../BatchJobs.R")

# load registry
reg <- BatchJobs::loadRegistry("modelEvalELRA-files/")

## get IDS to submit
# refit
id.fEt1t2C <- BatchExperiments::findExperiments(reg,
    prob.pattern = "dataFromExpertt1",
    algo.pattern = "gam",
    algo.pars    =(form %in%
        c("Expertt1", "Expertt2", "Statict1", "Statict2")))
BatchJobs::submitJobs(reg, ids=BatchJobs::findNotDone(reg, id.fEt1t2C))
BatchJobs::waitForJobs(reg)


## refit zero by t1, t2, L460
id.fZero <- BatchExperiments::findExperiments(reg,
    prob.pattern = "dataFromNoneNoneNone",
    algo.pattern = "gam",
    algo.pars    = (form %in% c("Expertt1", "Expertt2", "Statict1", "Statict2")))

BatchJobs::submitJobs(reg, ids=BatchJobs::findNotDone(reg, id.fZero))
BatchJobs::waitForJobs(reg)

# ## refit L460 with t1, t2, L460
id.fL460byt1t2L460 <- BatchExperiments::findExperiments(reg,
    prob.pattern = "dataFromStatict1",
    algo.pattern = "gam",
    algo.pars    = (form %in%
        c("Statict1", "Expertt1", "Expertt2", "Statict2")))

BatchJobs::submitJobs(reg, ids=sample(BatchJobs::findNotDone(reg, id.fL460byt1t2L460), 20))
BatchJobs::waitForJobs(reg)

## refit Expertt2 by Expertt1, Statict1, Expertt2, Statict2
id.Et2byEt1St1Et2 <- BatchExperiments::findExperiments(reg,
    prob.pattern = "dataFromExpertt2",
    algo.pattern = "gam",
    algo.pars    = (form %in%
        c("Statict1", "Expertt1", "Expertt2", "Statict2")))

BatchJobs::submitJobs(reg, ids=BatchJobs::findNotDone(reg, id.Et2byEt1St1Et2))
BatchJobs::waitForJobs(reg)
