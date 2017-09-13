#' Functions for simulation evaluation measures calculation
#'
#'
#' @rdname simulationA_measure
#' @keywords internal
#' @export
mse <- function(dif) {

    sum((dif)^2 )/length(dif)

}

#' @rdname simulationA_measure
#' @keywords internal
#' @export
mae <- function(dif) {

    sum(abs(dif))/length(dif)

}

#' @rdname simulationA_measure
#' @keywords internal
#' @export
getFreqOrigInSim <- function(df.fit) {

    bool.lo <- df.fit$fit.orig >= df.fit$lo
    bool.hi <- df.fit$fit.orig <= df.fit$hi
    bool.lo & bool.hi

}

#' @rdname simulationA_measure
#' @keywords internal
#' @export
getFreqSimInOrigKI <- function(df.fit) {

    in.lo.orig <- df.fit$fit >= df.fit$lo.orig
    in.hi.orig <- df.fit$fit <= df.fit$hi.orig

    in.lo.orig & in.hi.orig

}


#' extract comparison measures from data frame with comparisons differences from
#' simulations and true effects/true differences
#'
#' @rdname simulationA_measure
#' @import dplyr
#' @importFrom stats qnorm
#' @keywords internal
#' @export
get_stats <- function(comparisons.diff.df, ci.level=0.05) {

    clf <- tbl_df(comparisons.diff.df) %>%
        filter(!(fit==0 & lo==0 & hi==0)) %>%
        group_by(setting, form.int, subsetting, comparison, label)
    cd <- clf %>%
        summarise(
            RMSE     = round(sqrt(mean(fit.diff^2)), 3),
            MAE      = round(mean(abs(fit.diff)), 3),
            Coverage = round(mean(abs(fit.diff) <= qnorm(1-ci.level/2)*se), 3),
            belowCI  = round(mean(true.fit < fit - qnorm(1-ci.level/2)*se), 3),
            aboveCI  = round(mean(true.fit > fit + qnorm(1-ci.level/2)*se), 3))
    # return
    cd

}

#' @rdname simulationA_measure
#' @inherit get_stats
#' @keywords internal
#' @export
get_stats2 <- function(comparisons.diff.df, ci.level=0.05) {

    clf <- tbl_df(comparisons.diff.df) %>%
        filter(!(fit==0 & lo==0 & hi==0 & fit.diff==0)) %>%
        group_by(setting, form.int, subsetting, comparison, label, Lids)
    cd <- clf %>%
        summarise(
            RMSE     = round(sqrt(mean(fit.diff^2)), 3),
            MAE      = round(mean(abs(fit.diff)), 3),
            Coverage = round(mean(abs(fit.diff) <= qnorm(1-ci.level/2)*se), 3),
            belowCI  = round(mean(true.fit < fit - qnorm(1-ci.level/2)*se), 3),
            aboveCI  = round(mean(true.fit > fit + qnorm(1-ci.level/2)*se), 3))

    cd <- group_by(cd, setting, form.int, subsetting, comparison, label) %>%
        summarise(
            RMSE     = mean(RMSE),
            MAE      = mean(MAE),
            Coverage = mean(Coverage),
            belowCI  = mean(belowCI),
            aboveCI  = mean(aboveCI))
    # return
    cd

}

#' @rdname simulationA_measure
#' @keywords internal
#' @export
#' @import BatchJobs
#' @keywords internal
#' @export
get_valid_ids <- function(
    reg,
    prob.pattern,
    algo.pattern,
    algo.form,
    exclude.error   = TRUE,
    exclude.warning = FALSE) {

    id.form <- findExperiments(
        reg,
        prob.pattern = prob.pattern,
        algo.pattern = algo.pattern,
        algo.pars    = (form %in% algo.form))

    id.form  <- findDone(reg, ids=id.form)
    id.error <- getErrorMessages(reg, id.form)

    if(!all(is.na(id.error)) & exclude.error) {
        id.form <- setdiff(id.form, id.error)
    }
    id.warning <- grepLogs(reg, ids=id.form, pattern="warning")

    if(length(id.warning != 0) & exclude.warning) {
        id.form <- setdiff(id.form, id.warning)
    }

    # return
    id.form

}

#' @rdname simulationA_measure
#' @import BatchExperiments
#' @importFrom dplyr filter
#' @keywords internal
#' @export
get_experiments <- function(reg, probs, algos, forms, min.n) {

    experiments <- summarizeExperiments(
        reg,
        ids  = findDone(reg),
        show = c("prob", "algo", "form"))

    experiments <- filter(experiments, form %in% forms)
    if(!is.null(probs)) {
        experiments <- filter(experiments, prob %in% probs)
    }
    if(!is.null(algos)) {
        experiments <- filter(experiments, algo %in% algos)
    }
    experiments <- filter(experiments, .count >= min.n)

    # return
    experiments

}


#' Wrapper functions that obtains evaluated comparisons for jobs meeting
#' some criteria
#'
#' @rdname simulationA_measure
#' @import BatchJobs checkmate
#' @importFrom parallel mclapply
#' @importFrom dplyr left_join
#' @keywords internal
#' @export
get_all_comparisons <- function(
    reg,
    mc.cores        = 1,
    probs           = NULL,
    algos           = NULL,
    forms           = NULL,
    min.n           = 1,
    settings        = NULL,
    exclude.error   = TRUE,
    exclude.warning = TRUE) {

    experiments <- get_experiments(
        reg,
        probs = probs,
        algos = algos,
        forms = forms,
        min.n = min.n)


    stats.list <- mclapply(seq_len(nrow(experiments)),
        function(i) {

            prob.pattern <- experiments$prob[i]
            ids <- get_valid_ids(
                reg,
                prob.pattern    = prob.pattern,
                algo.pattern    = experiments$algo[i],
                algo.form       = experiments$form[i],
                exclude.error   = exclude.error,
                exclude.warning = exclude.warning)
            assert_numeric(ids, lower=1, finite=TRUE, any.missing=FALSE,
                min.len=1, unique=TRUE)
            # load results from registry
            lr <- lapply(ids, function(z) {
                    # print(z)
                    loadResult(reg, id=z)[["comparisons"]]
                })
            # ind.error <- sapply(lr, function(z) class(z)=="try-error")
            # lr <- lr[!ind.error]
            mcomp      <- melt_comparisons(lr)
            mcomp$prob <- factor(prob.pattern, levels=unique(experiments$prob))
            mcomp$algo <- factor(experiments$algo[i], levels=unique(experiments$algo))
            mcomp$form <- factor(experiments$form[i], levels=unique(experiments$form))
            # get true effects/comparisons
            truth.string <- sub("dataFrom", "", prob.pattern)
            # in case data was from null model (no nutrition effect), true
            # effect = 0
            if(grepl("NoneNoneNone", truth.string)) {
                truth <- data.frame(
                    intmid   = unique(mcomp$intmid),
                    true.fit = 0,
                    true.se  = 0,
                    true.hi  = 0,
                    true.lo  = 0)
            } else {
                truth <- get_true_diffs(
                    getProblem(reg, id=prob.pattern)$static,
                    truth.pattern=truth.string)
            }

            # add true effect and differences to truth to df
            mcomp <- add_truth(mcomp, truth)
            #add settings info
            if(!is.null(settings)) {
                dd <- left_join(mcomp, settings)
            } else {
                mcomp
            }

        }, mc.cores=mc.cores)

    return(stats.list)

}



#' @rdname simulationA_measure
#' @inherit get_all_comparisons
#' @keywords internal
#' @export
get_all_pvals <- function(
    reg,
    mc.cores        = 1,
    probs           = NULL,
    algos           = NULL,
    forms           = NULL,
    min.n           = 1,
    settings        = NULL,
    exclude.error   = TRUE,
    exclude.warning = TRUE) {

    experiments <- get_experiments(
        reg,
        probs = probs,
        algos = algos,
        forms = forms,
        min.n = min.n)

    pvals.list <- mclapply(seq_len(nrow(experiments)),
        function(i) {
            prob.pattern <- experiments$prob[i]

            ids <- get_valid_ids(
                reg,
                prob.pattern    = prob.pattern,
                algo.pattern    = experiments$algo[i],
                algo.form       = experiments$form[i],
                exclude.error   = exclude.error,
                exclude.warning = exclude.warning)

            assert_numeric(ids, lower=1, finite=TRUE, any.missing=FALSE,
                min.len=1, unique=TRUE)

            # load results from registry
            pvals <- lapply(ids, function(z) {
                rz <- unname(loadResult(reg, id=z)$summary.m$s.table[3:4, "p-value", drop=TRUE])
                data.frame(CII=rz[1], CIII=rz[2])
                })
            pvals.df      <- do.call(rbind, pvals)
            pvals.df$prob <- factor(prob.pattern, levels=unique(experiments$prob))
            pvals.df$algo <- factor(experiments$algo[i], levels=unique(experiments$algo))
            pvals.df$form <- factor(experiments$form[i], levels=unique(experiments$form))

            # get true effects/comparisons
            truth.string <- sub("dataFromF", "f", prob.pattern)
            pvals.df$truth <- truth.string

            # add settings info
            if(!is.null(settings)) {
                pvals.df <- dplyr::left_join(pvals.df, settings)
            }

            pvals.df

        }, mc.cores=mc.cores)


    return(pvals.list)

}
