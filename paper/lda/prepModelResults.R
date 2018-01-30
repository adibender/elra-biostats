library(tidyr)
library(dplyr)
library(ggplot2)
theme_set(theme_bw())
library(pammtools)
library(prodlim)
library(survival)
library(mgcv)
library(elrapack)

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


p1 <- ggplot_build(
  ggplot(baseline, aes(x = time, y = fit, ymin = lo, ymax = hi)) +
  geom_step() + geom_stepribbon(alpha = .2) + xlab("t") +
  ylab(expression(log(hat(lambda)[0](t)))) + ggtitle("Log-Baseline:"))
p2 <- ggplot_build(
  ggplot(apache, aes(x = time, y = fit, ymin = lo, ymax = hi)) +
  geom_step() + geom_stepribbon(alpha = .2) + xlab("t") +
  ylab(expression(hat(beta)[Apache](t))) +
  geom_hline(yintercept = 0, col="grey20", lty = 2) +
  ggtitle("Time-varying coefficient for Apache II:"))
p3 <- ggplot_build(
  ggplot(re, aes(sample = fit)) +
  geom_abline(slope = re$qqslope[1], intercept = re$qqintercept[1]) +
  geom_qq() + xlab("N(0,1) quantiles") +
  ylab(expression("ICU effects "~hat(gamma)[l[i]])) +
  ggtitle("QQ-Plot for ICU frailty:"))
p4 <- ggplot_build(
  ggplot(filter(pd, xlab == "Age"),
  aes(x = x, y = fit, ymin = low, ymax = high)) +
  geom_hline(yintercept = 0, col="grey20", lty = 2) +
  geom_line() + geom_ribbon(alpha = .2) +
  xlab("Age") + ylab(expression(hat(f)(Age))) +
  ggtitle("Varying coefficient of age:"))
p5 <- ggplot_build(
    ggplot(filter(pd, xlab == "BMI"),
  aes(x = x, y = fit, ymin = low, ymax = high)) +
  geom_hline(yintercept = 0, col="grey20", lty = 2) +
  geom_line() + geom_ribbon(alpha = .2) +
  xlab("BMI") +
  ylab(expression(hat(f)(BMI))) + ggtitle("Varying coefficient of BMI:"))

pdf("confoundert1.pdf", width=12, height=6)
cowplot::plot_grid(ggplot_gtable(p1), ggplot_gtable(p2), ggplot_gtable(p3),
  ggplot_gtable(p4), ggplot_gtable(p5), nrow=2)
dev.off()

################################################################################
## nutrition effects for Comparisons A-F (main analysis)

gg_maint1 <- ggplot(nutri.diffs, aes(x=intmid)) +
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

ggsave("maint1.pdf", gg_maint1, width = 9, height=6)
ggsave("maint1.eps", gg_maint1, width = 9, height=6)

## nutrition effects for Comparisons A-F (sensitivity analysis)
gg_maint1sens <- ggplot(nutri.diffs.sensitivity, aes(x=intmid)) +
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

ggsave("maint1sens.pdf", gg_maint1sens, width=9, height=6)
ggsave("maint1sens.eps", gg_maint1sens, width=9, height=6)

rm(mod.full.epxert.t1.sensitivity)
gc()

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

elra_heat <- ggplot(tidy) +
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
ggsave("ELRA_heatmaps.pdf", elra_heat, width = 10, height = 12)
ggsave("ELRA_heatmaps.eps", elra_heat, width = 10, height = 12)
ggsave("ELRA_heatmaps.jpg", elra_heat, width = 10, height = 12)

tidy_fit <- filter(tidy, estimate == "estimate")

ex1 <- data.frame(
  xmin = c(1:11)-0.5,
  xmax = c(2:12)-0.5,
  C    = rep(c("C2", "C3"), times = c(5, 6)),
  ymin = rep(-18, 11),
  ymax = rep(-19, 11))

elra_heat_est <- ggplot(tidy_fit) +
  geom_tile(aes(x=te, y=-t, fill = value, col=value)) +
  geom_contour(aes(x=te, y=-t, z =value), bins = 20, size=.5,
    colour = "gray70", ) +
  facet_wrap(~ C, nrow = 1) +
  coord_cartesian(xlim = c(1,11), ylim= -c(29.5, 4.5)) +
  scale_x_continuous(breaks = 1:11, minor_breaks = NULL,
    name = expression("Nutrition day "~t[e])) +
  scale_y_continuous(breaks = -(4:29 + .5), minor_breaks = NULL,
    name = expression(tilde(t)[j]),
    labels = paste0("(",4:29, ",", 5:30, "]"),
    sec.axis = dup_axis(name="")) +
  scale_colour_gradient2(na.value = 'grey90', guide = "none") +
  scale_fill_gradient2(na.value = 'grey90',
    name = expression("estimated "~tilde(h)(tilde(t)[j], t[e])),
    guide = guide_colourbar(direction = "vertical")) +
  theme(
    axis.text        = element_text(size = rel(1.4)),
    axis.title       = element_text(size = rel(1.5)),
    strip.text       = element_text(size = rel(1.5)),
    legend.text      = element_text(size = rel(1.4)),
    legend.title     = element_text(size = rel(1.5)),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    legend.position  = "right")

ggsave("elra_heat_estimation.pdf", elra_heat_est, width = 10, height = 6.5)
ggsave("elra_heat_estimation.jpg", elra_heat_est, width = 10, height = 6.5)

elra_heat_ex <- elra_heat_est+
  geom_rect(data = ex1, aes(ymin = ymin, ymax = ymax, xmax = xmax, xmin = xmin),
    alpha = 0, fill = "white", col = "black")

ggsave("elra_heat_ex.pdf", elra_heat_est, width = 10, height = 6.5)
ggsave("elra_heat_ex.jpg", elra_heat_est, width = 10, height = 6.5)


#### Figure point estimates with Lag-Lead window
#### combined figure with lag-lead window
tidy_est <- filter(tidy, estimate == "estimate")
gg_elra_est <- ggplot(tidy_est) +
  geom_tile(aes(x=te, y=-t, fill = value, col=value)) +
  geom_contour(aes(x=te, y=-t, z =value), bins = 20,
    colour = "gray70", lwd = .5) +
  facet_wrap(~ C, nrow = 1) +
  scale_colour_gradient2(na.value = 'grey90', guide = "none") +
  scale_fill_gradient2(na.value = 'grey90',
    name = expression("estimated "~tilde(h)(t, t[e]))) +
  theme(panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()) +
  # coord_cartesian(xlim = c(1,11), ylim= -c(29.5, 4.5)) +
  scale_x_continuous(breaks = 1:11, minor_breaks = NULL,
    name = expression("Nutrition day "~t[e])) +
  scale_y_continuous(breaks = -(4:29 + .5), minor_breaks = NULL,
    name = "t",
    labels = paste0("(",4:29, ",", 5:30, "]"),
    sec.axis = dup_axis(name=""))


# nutri        <- readRDS("../../dataCurrent/nutriOrigSmall.Rds")
seq.ints     <- seq_along(unique(nutri$intmid))
colindex     <- if(ncol(nutri$LHartlDynf)==12) -1 else 1:11
## graphic titles
dyn.title <- expression(
  atop("Dynamic time window",
    atop("lag = 4 days; lead = 4 + 2 x n days",
      "(n = number of days in the ICU)")))
static.title <- expression(atop("Static time window",
  atop("lag = 4 days; lead = 30 days")))


 ## hm: heat matrix containing the values to be ploted via heatmap()
int.names <- unique(elrapack::int_info(brks = c(0:30))$interval)
ll_dyn <- reshape2::melt(nutri$LHartlDynf[seq.ints, colindex])
ll_dyn$Var1 <- factor(ll_dyn$Var1, labels = int.names)
ll_dyn$Var2 <- factor(ll_dyn$Var2, labels = 1:11)
ll_dyn$type = "dynamic"

ll_stat <- reshape2::melt(nutri$LSimple460f[seq.ints, colindex])
ll_stat$Var1 <- factor(ll_stat$Var1, labels = int.names)
ll_stat$Var2 <- factor(ll_stat$Var2, labels = 1:11)
ll_stat$type = "static"

ll <- bind_rows(ll_dyn, ll_stat)

gg_dyn <- ggplot(ll_dyn, aes(x = Var2, y = rev(Var1))) +
    geom_tile(aes(fill = value), color = "lightgrey") +
    scale_fill_gradient(low = "whitesmoke", high = "grey20") +
    xlab("Protocol day") +
    scale_y_discrete("Interval j", labels = rev(int.names)) +
    theme(legend.position = "none") +
    ggtitle(dyn.title) +
    xlab(expression(paste("Nutrition day ", t[e])))

gg_stat <- ggplot(ll_stat, aes(x = Var2, y = rev(Var1))) +
    geom_tile(aes(fill = value), color = "lightgrey") +
    scale_fill_gradient(low = "whitesmoke", high = "grey20") +
    xlab("Protocol day") +
    scale_y_discrete("Interval j", labels = rev(int.names)) +
    theme(legend.position = "none") +
    ggtitle(static.title)  +
    xlab(expression(paste("Nutrition day ", t[e])))

library(cowplot)
ggsave("ELRA_heatmaps_LL.eps", plot_grid(gridExtra::grid.arrange(gg_dyn, gg_stat, nrow = 1), gg_elra_est, nrow = 2), width=9, height=12)
ggsave("ELRA_heatmaps_LL.pdf", plot_grid(gridExtra::grid.arrange(gg_dyn, gg_stat, nrow = 1), gg_elra_est, nrow = 2), width=9, height=12)
ggsave("ELRA_heatmaps_LL.jpg", plot_grid(gridExtra::grid.arrange(gg_dyn, gg_stat, nrow = 1), gg_elra_est, nrow = 2), width=9, height=12)

###################### model comparisons #######################################
#### compare via apparent C-Index
nutri <- readRDS("../../dataCurrent/nutri2.Rds")
mod_simple <- readRDS("../../runModelBatchJobs/simple_tdc_mod.Rds")
mod_simple2 <- readRDS("../../runModelBatchJobs/simple_tdc_modLLadequ.Rds")
mod_no_nutri <- readRDS("../../runModelBatchJobs/mod_no_nutrition.Rds")

## calculate apparent C-Index for all
cind_cumu            <- cindex_wrapper(mod.full.expert.t1, nutri)
cind_cumu$Model     <- "Main Model"
cind_simp            <- cindex_wrapper(mod_simple, nutri)
cind_simp$Model     <- "Model A1"
cind_simp2           <- cindex_wrapper(mod_simple2, nutri)
cind_simp2$Model    <- "Model A2"
cind_no_nutri        <- cindex_wrapper(mod_no_nutri, nutri)
cind_no_nutri$Model <- "No nutrition"

cind_mods <- bind_rows(cind_cumu, cind_simp, cind_simp2, cind_no_nutri)
saveRDS(cind_mods, "cindex_application.Rds")

gg_cind <- ggplot(cind_mods, aes(x = times, y = cindex, col = Model)) +
  geom_line() +
  ylim(c(0.5, 1)) +
  ylab("Apparent C-Index")

ggsave("cindexComparisonApplication.pdf", gg_cind, width = 8, height = 4)
ggsave("cindexComparisonApplication.jpg", gg_cind, width = 8, height = 4)

cind_mods %>%
  group_by(Model) %>%
  summarize(mcind = mean(cindex))



## objects containing information on nutrition history
protocols.df         <- readRDS("../../protocolsDF.Rds")
protocols.to.compare <- readRDS("../../protocolsToCompare.Rds")
med.patient          <- median_patient(subset(nutri, intmid==4.5))

## create "newdata" for all nutrtion profiles used in comparisons A - F
pat.list.full <- lapply(protocols.df, pattern_pat,
  ped = nutri, m = mod.full.expert.t1, median.patient = med.patient,
  effectname = "AdequacyCals")
pat.list.simple <- lapply(protocols.df, pattern_pat,
  ped = nutri, m = mod_simple, median.patient = med.patient,
  effectname = "adeq")
pat.list.simple2 <- lapply(protocols.df, pattern_pat,
  ped = nutri, m = mod_simple2, median.patient = med.patient,
  effectname = "LLadequ")
adeq_vars <- paste0("adeq", c("Tot0to30", "Tot30To70", "TotAbove70"))
pat.list.simple[[1]][, adeq_vars]
pat.list.simple[[2]][, adeq_vars]
adeq_vars2 <- paste0("LLadequ", c("Tot0to30", "Tot30To70", "TotAbove70"))
pat.list.simple2[[1]][, adeq_vars2]
pat.list.simple2[[4]][, adeq_vars2]

#### comparisons A - F
nutri_diffs_simple <- get_comp_diffs(
    pat.list.simple,
    protocols.df,
    protocols.to.compare,
    mod_simple,
    effectname = "adeq") %>%
  select(-L1, -se)  %>%
  mutate(Model = "Model A1")

nutri_diffs_simple2 <- get_comp_diffs(
    pat.list.simple2,
    protocols.df,
    protocols.to.compare,
    mod_simple2,
    effectname = "adeq") %>%
  select(-L1, -se)  %>%
  mutate(Model = "Model A2")

 nutri_diffs_full <- get_comp_diffs(
    pat.list.full,
    protocols.df,
    protocols.to.compare,
    mod.full.expert.t1,
    effectname = "AdequacyCals") %>%
  select(-L1, -se) %>%
  mutate(Model = "Main Model")

nutri_diffs_comp <- bind_rows(nutri_diffs_simple,
  nutri_diffs_simple2, nutri_diffs_full) %>%
  mutate(
    Model = factor(Model,
      levels = c("Main Model", "Model A1", "Model A2")))

gg_diffs_comp <- ggplot(nutri_diffs_comp, aes(x=intmid)) +
  geom_hline(yintercept = 0, lty=1, col="grey80") +
  geom_line(aes(y=fit, lty = Model, col = Model), lwd=1) +
  geom_ribbon(data = filter(nutri_diffs_comp, Model == "Main Model"),
    aes(ymin=lo, ymax=hi), alpha = 0.2) +
  facet_wrap(~label, nrow=2) +
  scale_y_continuous(
    breaks = log(sy <- c(0.25, 0.5, 0.75, 1, 1.25, 2, 4)),
    labels = sy,
    limits = c(-1.5, 1.5)) +
  scale_x_continuous(breaks = c(4, 10, 20, 30)) +
  ylab(expression(e[j])) +
  xlab(expression(t))+
  # Add Comparison info
  geom_text(
    dat = data.frame(
      x          = 10,
      y          = 1.3,
      label      = unique(nutri_diffs_simple$label),
      comparison = unique(nutri_diffs_simple$comparison)),
    aes(x=x, y=y, label=comparison), size=5) +
  theme(legend.position="bottom") +
  theme(
    axis.text    = element_text(size = rel(1.2)),
    axis.title   = element_text(size = rel(1)),
    title        = element_text(size = rel(1.3)),
    legend.title = element_text(size = rel(1)),
    legend.text  = element_text(size = rel(1.3)),
    strip.text   = element_text(size = rel(1.1)))

ggsave("main_comparison.jpg", gg_diffs_comp, width=9, height=6.5)
ggsave("main_comparison.pdf", gg_diffs_comp, width=9, height=6.5)

rm(mod_simple, mod_simple2, mod_no_nutri, mod.full.expert.t1)
gc()
