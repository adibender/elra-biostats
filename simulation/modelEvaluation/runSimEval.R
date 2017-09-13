library(BatchExperiments)
library(dplyr)
library(parallel)
library(checkmate)
# load registry and thereby necessary functions in:
# - evalMods.R
# - measureFunctions.R (get_stats)
reg <- BatchJobs::loadRegistry("../modelEvaluation/modelEvalELRA-files/")

experiments <- BatchExperiments::summarizeExperiments(
    reg,
    ids  = BatchJobs::findDone(reg),
    show = c("prob", "algo", "form"))

probs <- c("dataFromExpertt1", "dataFromExpertt2", "dataFromStatict1", "dataFromNoneNoneNone")
algos <- "gam"
forms <- c("Expertt1", "Expertt2", "Statict1", "Statict2")

experiments <- filter(experiments,
    form %in% forms,
    algo %in% algos,
    prob %in% probs)
experiments

n.prob <- length(unique(experiments$prob))
experiments$setting <- paste("Setting",
    c("I", "II", "IV", "III"))[rep(seq_len(n.prob), each=4)]
experiments <- arrange(experiments, setting)
experiments$setting <- as.factor(experiments$setting)
experiments$form.int <- factor(experiments$form,
    labels=letters[seq_along(unique(experiments$form))])
experiments <- transform(experiments,
    subsetting = interaction(setting, form.int))


# data frame containing new definition/ordering of comparisons
comp.new <- readRDS("../../protocolsToCompare.Rds")
save.path <- "../../paper/lda/simulation"

if(!test_directory_exists(save.path)) {
    dir.create(save.path)
}


## copparison stats: RMSE, Bias, Coverage
comp.list   <- get_all_comparisons(
    reg,
    mc.cores        = 16,
    min.n           = 1,
    probs           = c("dataFromExpertt1", "dataFromExpertt2",
        "dataFromStatict1", "dataFromNoneNoneNone"),
    algos           = c("gam"),
    forms           = c("Expertt1", "Expertt2", "Statict1", "Statict2"),
    settings        = experiments,
    exclude.warning = FALSE)

saveRDS(comp.list, "../../paper/lda/simulation/comparisonsList.Rds")


## pvalue stats (overall test)
pvals.standard <- get_all_pvals(
    reg             = reg,
    mc.cores        = 16,
    min.n           = 1,
    probs           = c("dataFromExpertt1", "dataFromExpertt2",
        "dataFromStatict1", "dataFromNoneNoneNone"),
    algos           = c("gam"),
    forms           = c("Expertt1", "Expertt2", "Statict1", "Statict2"),
    settings        = experiments,
    exclude.warning = FALSE)

saveRDS(pvals.standard, paste0(save.path, "/", "pvalsStandard.Rds"))
