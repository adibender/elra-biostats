#' Functions for evaluation of Simulation Study Part A.
#'
#' The main function \code{eval_model} is called in the wrapper algorithm after
#' model is fit. Other functions basically helper functions.
#'
#' @rdname eval_model
#' @param pam A piece-wise exponential model (gam).
#' @param X.list List of model matirces representing protocol Comparisons.
#' @param protocols.to.compare  Data frame of protocols that should be compared
#' @param protocols.df Dayly protocol regimen for protocols defined in \code{protocols.to.compare}.
#' @param maxdays.tdc Maximal number of days nutrition protocol was recorded.
#' @param intmid Interval midpoints of intervals in which follow-up is devided in.
#' @param ... further arguments, currently not used.
#' @importFrom stats coefficients
#'
#' @return  List containing mdoel coefficients and their covariance, summary of
#' the model object, Estimated ln HR of the comparisons and p-values for each
#' comparison.
#'
#' @keywords internal
#' @export
eval_model <- function(
    pam,
    X.list,
    protocols.to.compare,
    protocols.df,
    maxdays.tdc,
    intmid,
    debug = TRUE, ...) {

    pam.coefs          <- coefficients(pam)
    pam.Vp             <- pam$Vp
    pam.summary        <- summary(pam)
    pam.summary$family <- NULL # otherwise uses up way more space on saving
    attr(pam.summary, ".Environment") <- NULL # dito

    comparisons <- get_comparison_diffs(
        coef.vec             = pam.coefs,
        Vp                   = pam.Vp,
        X.list               = X.list,
        protocols.to.compare = protocols.to.compare,
        protocols.df         = protocols.df,
        maxdays.tdc          = maxdays.tdc,
        intmid               = intmid)

    pvals.comparisons <- get_comparison_pvals(
        model                = pam,
        X.list               = X.list,
        protocols.df         = protocols.df,
        protocols.to.compare = protocols.to.compare,
        maxdays.tdc          = maxdays.tdc,
        debug                = debug)

    # return
    list(
        coefs       = pam.coefs,
        Vp          = pam.Vp,
        comparisons = comparisons,
        convergence = pam$converged,
        edf         = pam$edf,
        edf1        = pam$edf1,
        summary.m   = pam.summary,
        p.vals      = pvals.comparisons)
}

#' @rdname eval_model
#' @importFrom reshape2 melt
#' @keywords internal
#' @export
get_comparison_diffs <- function(
    coef.vec,
    Vp,
    X.list,
    protocols.to.compare,
    protocols.df,
    maxdays.tdc,
    intmid) {

    diff.labs <- apply(protocols.to.compare, 1,
        function(z) {
            paste(
                pattern_label(protocols.df[, z[1]], maxdays.tdc = maxdays.tdc),
                "vs. \n",
                pattern_label(protocols.df[, z[2]], maxdays.tdc = maxdays.tdc))
        })
    effect.diffs <- apply(protocols.to.compare, 1,
        function(z) {
            diff.i <- as.data.frame(
                get_effect_diffs(
                    coef.vec,
                    Vp,
                    X1     = X.list[[z[1]]],
                    X2     = X.list[[z[2]]],
                    intmid = intmid))
        })

    names(effect.diffs) <- diff.labs
    m.effect.diffs <- melt(
        effect.diffs,
        id.vars=c("fit", "intmid", "se", "hi", "lo"))
    m.effect.diffs$label <- factor(m.effect.diffs$L1, levels = diff.labs)
    m.effect.diffs$L1    <- NULL

    ## return
    m.effect.diffs

}

#' @rdname eval_model
#' @keywords internal
#' @export
get_df <- function(model, useColumns) {

    df <- sum(model$edf[useColumns])
    if (!is.null(model$edf1)) df <- sum(model$edf1[useColumns])
    df <- min(length(useColumns), df)

    # return
    df
}

#' @rdname eval_model
#' @keywords internal
#' @export
get_custom_df <- function(edf, X1, X2) {
    edf2 <- X2 %*% t(t(edf))
    edf1 <- X1 %*% t(t(edf))

}

#' @rdname eval_model
#' @keywords internal
#' @export
get_rdf <- function(model) {

    if (model$scale.estimated) {
        rdf <-  length(model$y) - sum(model$edf)
    } else rdf <- -1

    # return
    rdf
}

#' @rdname eval_model
#' @keywords internal
#' @export
test_effect_diffs <- function(
    model,
    X1,
    X2,
    effectname = "AdequacyCals",
    debug      = FALSE,
    debugX     = FALSE) {

    # extract infos
    coef.vec <- coefficients(model)
    Vp <- model$Vp

    # sometimes coef vec shorter than dim of covarianz matrix
    useColumns <- grep(effectname, names(coef.vec))
    if(length(coef.vec)!=ncol(X1)) {
        useColumns.X <- grep(effectname, colnames(X1))
    } else {
        useColumns.X <- useColumns
    }


    covCoefs   <- Vp[useColumns, useColumns]
    X1         <- X1[,useColumns.X]
    X2         <- X2[,useColumns.X]
    X          <- X2 - X1

    if(debugX) return(list(X1=X1, X2=X2))

    dropcols   <- which(apply(X, 2, function(x) all(x == 0)))

    if(length(dropcols) > 0) {
        useColumns <- useColumns[-dropcols]
        X          <- X[,-dropcols]
        covCoefs   <- covCoefs[-dropcols, -dropcols]
    }

    df  <- get_df(model, useColumns)
    rdf <- get_rdf(model)

    # return
    if(debug) {
        return(list(
            p      = coef.vec[useColumns],
            X      = X,
            V      = covCoefs,
            rank   = df,
            type   = 0,
            res.df = rdf
            ))
    } else {
        mgcv:::testStat(
            p      = coef.vec[useColumns],
            X      = X,
            V      = covCoefs,
            rank   = df,
            type   = 0,
            res.df = rdf)
    }

}

#' @rdname eval_model
#' @keywords internal
#' @export
get_comparison_pvals <- function(
    model,
    X.list,
    protocols.df,
    protocols.to.compare,
    maxdays.tdc,
    debug=FALSE) {

    comp.pval <- apply(protocols.to.compare, 1,
        function(z) {
            cat(z[1], "vs.", z[2], "\n")
            test_effect_diffs(
                model = model,
                X1    = X.list[[z[1]]],
                X2    = X.list[[z[2]]],
                debug = debug)
        })
    if(debug) {
        return(comp.pval)
    }
    else {
        comp.pval <- comp.pval$pval
    }

    labs <- apply(protocols.to.compare, 1,
        function(z) {
            paste(
                pattern_label(protocols.df[, z[1]], maxdays.tdc),
                "vs. \n",
                pattern_label(protocols.df[, z[2]], maxdays.tdc))
        })

    pval.df <- data.frame(
        label  = labs,
        pval   = round(comp.pval, 4))
    pval.df$pval.num <- pval.df$pval
    pval.df$pval <- with(pval.df,
        ifelse(pval > 0.001, paste0("p = ", pval), "p < 0.001"))

    ## return
    pval.df

}

#' @rdname eval_model
#' @param effectname Name of variable(s) for which the effect should be extracted.
#' \code{"AdequacyCals"} or \code{"ProtCat"}.
#' @param X1 Design matrix.
#' @param X2 Design matrix.
#' @importFrom stats qnorm
#' @keywords internal
#' @export
get_effect_diffs <- function(
    coef.vec,
    Vp,
    X1,
    X2,
    intmid,
    effectname = "AdequacyCals",
    alpha      = 0.05,
    debug      = FALSE) {

    # sometimes coef vec shorter than dim of X (depending on simulation setting)
    # X has always dimension for max/full model
    useColumns <- grep(effectname, names(coef.vec))
    if(length(coef.vec)!=ncol(X1)) {
        useColumns.X <- grep(effectname, colnames(X1))
    } else {
        useColumns.X <- useColumns
    }

    covCoefs   <- Vp[useColumns, useColumns]
    X1         <- X1[,useColumns.X]
    X2         <- X2[,useColumns.X]
    X          <- X2 - X1
    if(debug) return(X)
    fit        <- drop(X %*% coef.vec[useColumns])
    se         <- sqrt(rowSums((X %*% covCoefs) * X))
    hi         <- fit + qnorm(1-alpha/2) * se
    lo         <- fit - qnorm(1-alpha/2) * se

    #return
    cbind.data.frame(
        intmid = intmid,
        fit    = fit,
        se     = se,
        hi     = hi,
        lo     = lo)

}

#' @rdname eval_model
#' @importFrom stats coef
#' @keywords internal
#' @export
proto_lp <- function(
    pattern,
    model,
    ped,
    median.x,
    effectname = "AdequacyCals",
    time.var   = "intmid",
    id.var     = "CombinedID",
    type       = "lpmatrix",
    se.fit     = FALSE) {

    # pick a patient that remained under risk all the way to overwrite their data:
    proto.x <- {
        ind <- which(ped[[time.var]] == max(ped[[time.var]]))[1]
        id <- ped[[id.var]][ind]
        subset(ped, ped[[id.var]] == id)
    }
    # overwrite with median data for confounders:
    for(var in colnames(median.x)) {
        proto.x[, var] <- median.x[1, var]
    }

    # create TDC covariates, given TDC pattern
    gc <- grep("AdequacyCals", names(coef(model)), value = TRUE)

    if(length(gc) == 0) {
        # no usefull comparison, if no TDC coefs in model
        return(NULL)
    } else {
        effectname <- sub(".*I\\(", "", gc)
        effectname <- unique(sub("Tot.*", "", effectname))

        low.var    <- paste0(effectname, "Tot0to30")
        mid.var    <- paste0(effectname, "Tot30To70")
        full.var   <- paste0(effectname, "TotAbove70")
        # switch
        proto.x[[low.var]]  <- outer(rep(1, nrow(proto.x)), pattern == "low")
        proto.x[[mid.var]]  <- outer(rep(1, nrow(proto.x)), pattern == "mid")
        proto.x[[full.var]] <- outer(rep(1, nrow(proto.x)), pattern == "full")

        # TDC covariates are wrongly created if any entry != 1
        TDC.all <- proto.x[[low.var]] + proto.x[[mid.var]] + proto.x[[full.var]]
        if(!all(TDC.all==1)) {
            stop("Error in creation of TDC covariates")
        }

    }

    # correct lag-lead (and everything else) is automatically selected b/c model
    # objec is provided
    predict(model, newdata=proto.x, type=type, se.fit=se.fit)

}

#' @rdname eval_model
#' @keywords internal
#' @export
get_true_diffs <- function(
    static,
    truth.pattern = "Expertt1",
    rename.cols   = c("fit", "se", "hi", "lo")) {

    true.diffs <- get_comparison_diffs(
        static$coefs.list[[truth.pattern]],
        static$Vp.list[[truth.pattern]],
        get_X(static$Xcompare.list, truth.pattern),
        static$protocols.to.compare,
        static$protocols.df,
        static$maxdays.tdc,
        static$intmid)

    rename.ind <- colnames(true.diffs) %in% rename.cols
    colnames(true.diffs)[rename.ind] <- paste("true", rename.cols, sep=".")

    # return
    true.diffs

}

#' @rdname eval_model
#' @keywords internal
#' @export
get_X <- function(Xcompare.list, truth.pattern="Expertt1") {

    ## Needed data set is the same for most sensitivity analyses, therefore
    # Xcompare.list[1] can be used
    if(!(truth.pattern %in% names(Xcompare.list)))
        truth.pattern <- names(Xcompare.list)[1]

    Xcompare.list[[truth.pattern]]

}

#' @rdname eval_model
#' @importFrom reshape2 melt
#' @keywords internal
#' @export
melt_comparisons <- function(
    comparisons.list,
    id.vars=c("fit", "intmid", "se", "hi", "lo", "label")) {

    # return
    mcomp <- melt(comparisons.list, level="ids", id.vars=id.vars)
    mcomp$Lids <- factor(mcomp$Lids,
        labels=paste0("JobID", unique(mcomp$Lids)))
    mcomp$comparison <- factor(mcomp$label,
        label=paste0("Comparison ", LETTERS[seq_len(nlevels(mcomp$label))]))

    return(mcomp)

}


#' @rdname eval_model
#' @importFrom dplyr left_join
#' @keywords internal
#' @export
add_truth <- function(comparisons.df, true.comparisons) {

    comparisons.df <- left_join(comparisons.df, true.comparisons)
    comparisons.df$fit.diff <- comparisons.df$fit - comparisons.df$true.fit

    # return
    comparisons.df

}