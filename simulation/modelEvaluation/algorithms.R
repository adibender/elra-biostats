##' will be added as "algorithm" to a BatchExperiments registry
##' @param static Is loaded at the beginning and provided to all "problems" and
##' "algorithms" that will be added to the registry
##' @param dynamic A simulated data set. Output of the function provided as
# \code{dynamic} argument to \code{addProblem}.
##' @param form Formula object, that will be passed to \code{bam}, which then
# fits a model on \code{dynamic}.


bam_wrapper <- function(static, dynamic, form, maxit=20, ...) {

    # fit model on current data set (dynamic)
    mod <- bam(
        formula = static[["formulas.list"]][[form]],
        data    = dynamic,
        family  = poisson(),
        offset  = offset,
        control = gam.control(trace = TRUE, maxit = maxit), ...)

    # obtain correct list of comparisons for given form(ula)
    X.list <- get_X(static$Xcompare.list, form)
    # call wrapper to obtain summary statistics for the current model
    # needed later for plots and to obtain comparisons to true model
   eval_model(
        pam                  = mod,
        X.list               = X.list,
        protocols.df         = static$protocols.df,
        protocols.to.compare = static$protocols.to.compare,
        maxdays.tdc          = static$maxdays.tdc,
        intmid               = static$intmid)

}

gam_wrapper <- function(static, dynamic, form, debug=TRUE, maxit=20, ...) {

    # fit model on current data set (dynamic)
    mod <- gam(
        formula = static[["formulas.list"]][[form]],
        data    = dynamic,
        family  = poisson(),
        offset  = offset,
        method  = "REML",
        control = gam.control(trace = TRUE, maxit = maxit), ...)

    # obtain correct list of comparisons for given form(ula)
    X.list <- get_X(static$Xcompare.list, form)
    # call wrapper to obtain summary statistics for the current model
    # needed later for plots and to obtain comparisons to true model
   eval_model(
        pam                  = mod,
        X.list               = X.list,
        protocols.df         = static$protocols.df,
        protocols.to.compare = static$protocols.to.compare,
        maxdays.tdc          = static$maxdays.tdc,
        intmid               = static$intmid,
        debug                = debug)

}
