library(tidyr)
library(dplyr)
library(ggplot2)
theme_set(theme_bw())

model.path      <- "../../runModelBatchJobs/gamBatch/results/"
model.path.hosp <- "../../runModelBatchJobs/gamHospBatch/results/"
mod.name        <- "fullExpertt1Calories.Rds"

## load main model and sensitivity model
mod.full.expert.t1             <- readRDS(paste0(model.path, mod.name))
mod.full.epxert.t1.sensitivity <- readRDS(paste0(model.path.hosp, mod.name))

## load data sets for main and sensitivity analysis
nutri                <- readRDS("../../dataCurrent/nutriOrigSmall.Rds")
nutri.hosp           <- readRDS("../../dataCurrentHosp/nutriOrigSmall.Rds")
## objects containing information on nutrition history
protocols.df         <- readRDS("../../protocolsDF.Rds")
protocols.to.compare <- readRDS("../../protocolsToCompare.Rds")


## data frame with one row containing covariate information for "median patient"
med.patient      <- median_patient(subset(nutri, intmid==4.5))
med.patient.hosp <- median_patient(subset(nutri.hosp, intmid==4.5))

pat.list <- lapply(protocols.df, pattern_pat,
  ped = nutri, m = mod.full.expert.t1, median.patient = med.patient)
pat.list.hosp <- lapply(protocols.df, pattern_pat,
  ped = nutri.hosp, m = mod.full.epxert.t1.sensitivity, median.patient = med.patient.hosp)

## calculate log-hazard differences for custom comparisons A-F
nutri.diffs <- get_comp_diffs(pat.list, protocols.df, protocols.to.compare,
  mod.full.expert.t1) %>% select(-L1, -se) %>%
	gather(type, value, hi, lo)

nutri.diffs.sensitivity <- get_comp_diffs(pat.list.hosp, protocols.df, protocols.to.compare,
  mod.full.epxert.t1.sensitivity) %>% select(-L1, -se) %>%
	gather(type, value, hi, lo)

################################################################################
## visualize model effects
library(pam)
ints <- sort(unique(mod.full.expert.t1$model$intmid))
pd <- mod.full.expert.t1$model[rep(1, length(ints)),]
pd$ApacheIIScore <- 1
pd$intmid <- ints
baseline_data <- predict(mod.full.expert.t1, pd, type = "terms", se = TRUE)
baseline <- data_frame(
  t = ints,
  fit = baseline_data$fit[, "s(intmid)"] + coef(mod.full.expert.t1)[1],
  se = baseline_data$se.fit[, "s(intmid)"] + sqrt(vcov(mod.full.expert.t1)[1,1])) %>%
  mutate(tstart = floor(t) + 0.01, tstop = ceiling(t) - 0.01,
    hi = fit + 2*se, lo = fit - 2*se) %>%
  gather(key, value = "time", tstart, tstop)
apache_cols <- grep("ApacheIIScore", colnames(baseline_data$fit))
apache <- data_frame(
  t = ints,
  fit = rowSums(baseline_data$fit[, apache_cols]),
  se = rowSums(baseline_data$se.fit[, apache_cols])) %>%
  mutate(tstart = floor(t) + 0.01, tstop = ceiling(t) - 0.01,
    hi = fit + 2*se, lo = fit - 2*se) %>%
  gather(key, value = "time", tstart, tstop)

pd <- tidy_smooth(mod.full.expert.t1)
re <- tidy_re(mod.full.expert.t1)

p1 <- ggplot_gtable(ggplot_build(
  ggplot(baseline, aes(x = time, y = fit, ymin = lo, ymax = hi)) +
  geom_step() + geom_stepribbon(alpha = .2) + xlab("t") +
  ylab(expression(log(hat(lambda)[0](t)))) + ggtitle("Log-Baseline:")))
p2 <- ggplot_gtable(ggplot_build(
  ggplot(apache, aes(x = time, y = fit, ymin = lo, ymax = hi)) +
  geom_step() + geom_stepribbon(alpha = .2) + xlab("t") +
  ylab(expression(hat(beta)[Apache](t))) +
  geom_hline(yintercept = 0, col="grey20", lty = 2) +
  ggtitle("Time-varying coefficient for Apache II:")))
p3 <- ggplot_gtable(ggplot_build(
  ggplot(re, aes(sample = fit)) +
  geom_abline(slope = re$qqslope[1], intercept = re$qqintercept[1]) +
  geom_qq() + xlab("N(0,1) quantiles") +
  ylab(expression("ICU effects "~hat(gamma)[l[i]])) +
  ggtitle("QQ-Plot for ICU frailty:")))
p4 <- ggplot_gtable(ggplot_build(
  ggplot(filter(pd, xlab == "Age"),
  aes(x = x, y = fit, ymin = low, ymax = high)) +
  geom_hline(yintercept = 0, col="grey20", lty = 2) +
  geom_line() + geom_ribbon(alpha = .2) +
  xlab("Age") + ylab(expression(hat(f)(Age))) +
  ggtitle("Varying coefficient of age:")))
p5 <- ggplot_gtable(ggplot_build(
    ggplot(filter(pd, xlab == "BMI"),
  aes(x = x, y = fit, ymin = low, ymax = high)) +
  geom_hline(yintercept = 0, col="grey20", lty = 2) +
  geom_line() + geom_ribbon(alpha = .2) +
  xlab("BMI") +
  ylab(expression(hat(f)(BMI))) + ggtitle("Varying coefficient of BMI:")))

pdf("confoundert1.pdf", width=12, height=6)
cowplot::plot_grid(p1, p2, p3, p4, p5, nrow=2)
dev.off()

################################################################################
## nutrition effects for Comparisons A-F (main analysis)
pdf("maint1.pdf", width=9, height=6)
ggplot(nutri.diffs, aes(x=intmid)) +
  geom_hline(yintercept = 0, lty=1, col="grey70") +
	geom_line(aes(y=fit), lwd=1) +
	geom_line(aes(y=value, group=type), lty=2, lwd=1) +
	facet_wrap(~label, nrow=2) +
  scale_y_continuous(
    breaks = log(sy <- c(0.25, 0.5, 0.75, 1, 1.25, 2, 4)),
    labels = sy,
    limits = c(-1.5, 1.5)) +
  ylab(expression(e[j])) +
  xlab(expression(t))+
  # Add Comparison info
  geom_text(
  	dat = data.frame(
      x          = 10,
      y          = 1.3,
      label      = unique(nutri.diffs$label),
      comparison = unique(nutri.diffs$comparison)),
  	aes(x=x, y=y, label=comparison), size=5) +
  theme(legend.position="none") +
  theme(
    axis.text    = element_text(size = rel(1.2)),
    axis.title   = element_text(size = rel(1)),
    title        = element_text(size = rel(1.3)),
    legend.title = element_text(size = rel(1)),
    legend.text  = element_text(size = rel(1.3)),
    strip.text   = element_text(size = rel(1.1)))
dev.off()

## nutrition effects for Comparisons A-F (sensitivity analysis)
pdf("maint1sens.pdf", width=9, height=6)
ggplot(nutri.diffs.sensitivity, aes(x=intmid)) +
  geom_hline(yintercept = 0, lty=1, col="grey70") +
	geom_line(aes(y=fit), lwd=1) +
	geom_line(aes(y=value, group=type), lty=2, lwd=1) +
	facet_wrap(~label, nrow=2) +
  scale_y_continuous(
    breaks = log(sy <- c(0.25, 0.5, 0.75, 1, 1.25, 2, 4)),
    labels = sy,
    limits = c(-1.5, 1.5)) +
  geom_hline(yintercept=0, lty=2) +
  ylab(expression(e[j])) +
  xlab(expression(t)) +
  # Add Comparison info
  geom_text(
  	data=data.frame(
      x          = 10,
      y          = 1.3,
      label      = unique(nutri.diffs.sensitivity$label),
      comparison = unique(nutri.diffs.sensitivity$comparison)),
  	aes(x=x, y=y, label=comparison), size=5) +
  theme(legend.position="none") +
  theme(
    axis.text    = element_text(size = rel(1.2)),
    axis.title   = element_text(size = rel(1)),
    title        = element_text(size = rel(1.3)),
    legend.title = element_text(size = rel(1)),
    legend.text  = element_text(size = rel(1.3)),
    strip.text   = element_text(size = rel(1.1)))
dev.off()

################################################################################
## heat map for partial nutrition effects

pd <- mod.full.expert.t1$model[1,]
ints <- sort(unique(mod.full.expert.t1$model$intmid))
days <- 1:11
pd <- pd[rep(1, length(ints)*length(days)),]
pd$DaysMat <- rep(days, times = length(ints))
pd$intMat <- rep(ints, each = length(days))
pd[["I(AdequacyCalsTotAbove70 * LHartlDynf)"]] <-
  pd[["I(AdequacyCalsTot30To70 * LHartlDynf)"]] <- 1
preds <- predict(mod.full.expert.t1, newdata = pd, type = "terms", se = TRUE)
use <- grep("Adequacy", colnames(preds$fit))
fit <- preds$fit[, use]; colnames(fit) <-  c("C_III", "C_II")
hi <- fit + 2 * preds$se.fit[, use]; colnames(hi) <- paste0("hi_", colnames(fit))
lo <- fit - 2 * preds$se.fit[, use]; colnames(lo) <- paste0("lo_", colnames(fit))
window <- function(t, te) ifelse(t > te + 4 &
    floor(t) < te + 3 + 2*(te + 4), 1, NA)
tidy <- cbind(select(pd, intMat, DaysMat) %>% rename(t = intMat, te = DaysMat),
    fit, hi, lo) %>% gather(key, value, -t, -te)  %>%
  mutate(value = window(t, te) * value,
    key = factor(key,
      labels = c("es_C2", "es_C3", "hi_C2", "hi_C3", "lo_C2", "lo_C3"))) %>%
  separate(key, into = c("estimate", "C"), sep = 3) %>%
  mutate(estimate =
      factor(estimate, labels = c("estimate", "upper CI", "lower CI")))

pdf("ELRA_heatmaps.pdf", width=10, height=12)
ggplot(tidy) +
  geom_tile(aes(x=te, y=-t, fill = value, col=value)) +
  geom_contour(aes(x=te, y=-t, z =value), bins = 20, size=.1,
    colour = "gray70") +
  facet_grid(estimate ~ C) +
  coord_cartesian(xlim = c(1,11), ylim= -c(29.5, 4.5)) +
  scale_x_continuous(breaks = 1:11, minor_breaks = NULL,
    name = expression("Nutrition day "~t[e])) +
  scale_y_continuous(breaks = -(4:29 + .5), minor_breaks = NULL,
    name = "t",
    labels = paste0("(",4:29, ",", 5:30, "]"),
    sec.axis = dup_axis(name="")) +
  scale_colour_gradient2(na.value = 'grey90', guide = "none") +
  scale_fill_gradient2(na.value = 'grey90',
    name = expression("estimated "~tilde(h)(t, t[e])),
    guide = guide_colourbar(direction = "horizontal")) +
  theme(panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    legend.position = "bottom")
dev.off()