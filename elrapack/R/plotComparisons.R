#' Plot functions for hazard ratio comparisons given different nutrition profiles
#'
#' \code{ggplot_comparison} is the base function. Other functions mostly extensions.
#'
#' @rdname comparisons
#'
#' @import ggplot2
#' @importFrom dplyr select
#' @param comparisons A data frame defining the comparisons of interest.
#' @param true.comparisons If provided, these true comparisons will be added to
#' the plot of simulated comparisons
#' @keywords internal
#' @export
ggplot_comparison <- function(comparisons, true.comparisons=NULL, ...) {

    if(!is.data.frame(comparisons)) {
        comparisons <- melt_comparisons(comparisons)
    }
    if(!is.null(true.comparisons)) {
        comparisons <- add_truth(comparisons, true.comparisons)
    }

    comp <- unique(dplyr::select(comparisons, label, comparison))
    comp$intmid <- 17.5
    comp$fit <- -1.5
    ggplot(comparisons, aes(x = intmid, y = fit)) +
        # basics and labels
        ylab("log-hazard difference") +
        xlab("Days after ICU admission") +
        theme(legend.position = "none") +
        # scale_x_continuous(minor_breaks = c(4, int_info()$tend)) +
        theme(
            axis.text    = element_text(size = rel(1.2)),
            axis.title   = element_text(size = rel(1)),
            title        = element_text(size = rel(1.3)),
            legend.title = element_text(size = rel(1)),
            legend.text  = element_text(size = rel(1.3)),
            strip.text   = element_text(size = rel(1.1))) +
        ylim(c(-2.25, 2.25)) +
        facet_wrap(~label, nrow=2) +
        geom_text(data=comp, hjust=1, vjust=0.5,
            aes(x=intmid, y=fit, label=comparison))

}

#' @rdname comparisons
#' @inherit ggplot_comparison
#' @keywords internal
#' @export
ggbox_comparisons <- function(gg.base.comparison, ...) {

    gg.base.comparison + #aes(factor(intmid)) +
        geom_hline(yintercept=0, lty=2, col=1) +
        geom_boxplot(aes(group=intmid), outlier.size=1) +
        geom_line(aes(y=true.fit))

}

#' @rdname comparisons
#' @inherit ggplot_comparison
#' @keywords internal
#' @export
gglines_comparisons <- function(gg.base.comparison, alpha=0.1, ...) {

    gg.base.comparison +
        geom_line(aes(group=Lids), col="grey70", alpha=alpha) +
        geom_line(aes(y=true.fit), col=1, lwd=1.2)

}

#' @rdname comparisons
#' @inherit ggplot_comparison
#' @keywords internal
#' @export
ggbox_comparison_diffs <- function(gg.base.comparison, ...) {

    gg.base.comparison +
        aes(y=fit.diff) +
        ylab("simulated - true") +
        geom_boxplot(aes(group=intmid), outlier.size=1) +
        geom_hline(yintercept = 0, lty=2)

}

#' @rdname comparisons
#' @inherit ggplot_comparison
#' @keywords internal
#' @export
gglines_comparison_diffs <- function(gg.base.comparison, alpha=0.1, ...) {

    gg.base.comparison +
        aes(y=fit.diff) +
        ylab("simulated - true") +
        geom_line(aes(group=Lids), alpha=alpha, col="grey70") +
        geom_hline(yintercept = 0, lty=2)

}

#' @rdname comparisons
#' @inherit ggplot_comparison
#' @keywords internal
#' @export
ggtext_simstats <- function(simstats, id.vars=c("label", "comparison"), ...) {

    mstats <- melt(simstats, id.vars=id.vars)
    mstats <- filter(mstats, variable %in% c("RMSE", "Coverage"))
    mstats$intmid <- NA
    mstats$intmid[mstats$variable=="RMSE"] <- 30
    # mstats$intmid[mstats$variable=="MAE"] <- 60
    mstats$intmid[mstats$variable=="Coverage"] <- 30
    mstats$fit[mstats$variable=="RMSE"] <- 1.2
    # mstats$fit[mstats$variable=="MAE"] <- 1.5
    mstats$fit[mstats$variable=="Coverage"] <- 0.85
    mstats$stats <- paste0(mstats$variable, "=", round(mstats$value, 3))

    # return
    geom_text(data=mstats, hjust=1, size=5, vjust=0,
        aes(x=intmid, y=fit, group=variable, label=stats))

}

#' @rdname comparisons
#' @inherit ggplot_comparison
#' @keywords internal
#' @export
ggcomparison_lines2 <- function(setting, alpha=0.1) {

    gg.base <- ggplot_comparison(setting)
    simstats <- get_stats2(gg.base$data)
    simstats$form.int <- NULL
        # create text layer with summary statistics
    if(!is.null(simstats[["setting"]])) {
        id.vars <- c("label", "comparison", "setting", "subsetting")
    } else {
        id.vars <- c("label", "comparison")
    }
    gg.simstats <- ggtext_simstats(simstats, id.vars=id.vars)
    gg.base     <- gg.base + gg.simstats

    gglines_comparisons(gg.base, alpha=alpha) +
    scale_y_continuous(
        breaks = log(sy <- c(0.25, 0.5, 0.75, 1, 1.25, 2, 4)),
        labels = sy) +
    coord_cartesian(ylim=c(-1.5, 1.5)) +
    geom_hline(yintercept=0, lty=2) +
    ylab(expression(e[j])) +
    xlab(expression(t))

}