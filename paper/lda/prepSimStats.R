## libraries
library(tables) # tabular function for creation of complex latex formated tables
library(grid) # combine multiple ggplots
library(gridExtra) # dito
library(ggplot2) # used to create graphics
theme_set(theme_bw()) # sets default ggplot graphics theme to black/white
library(BatchExperiments) # contains functions to extract information from registry containing simulation results
library(reshape2) # transform data sets to ggplot2 exploitable format
library(dplyr) # general purpose data manipulation

## data sets
# comparisonsList.Rds is created in "../../simulation/modelEvaluation/runSimEval.R"
# contains cumulative effects as log-hazard-differences (of 6 Comparisons: A-F)
cl <- readRDS("simulation/comparisonsList.Rds")

## load the registry to obtain information regarding different simulation settings
reg <- BatchJobs::loadRegistry("../../simulation/modelEvaluation/modelEvalELRA-files/")

# experiments/simulations are defined by data simulation ("prob"), evaluation
# algorithm ("prob", here this is always the same, namely mgcv::gam), and
# the model specification ("form")
experiments <- BatchExperiments::summarizeExperiments(
  reg,
  ids  = BatchJobs::findDone(reg),
  show = c("prob", "algo", "form"))

# extract experiments/simulations that will be evaluated (here all)
probs <- c("dataFromExpertt1", "dataFromExpertt2", "dataFromStatict1", "dataFromNoneNoneNone")
algos <- "gam"
forms <- c("Expertt1", "Expertt2", "Statict1", "Statict2")
experiments <- filter(experiments,
  form %in% forms,
  algo %in% algos,
  prob %in% probs)

## Add some nicer labeling and Setting information
# Later used as labels for  axes/legends/panes
n.prob              <- length(unique(experiments$prob))
experiments$setting <- paste("Setting", c("I", "II", "IV", "III"))[rep(seq_len(n.prob),
  each=4)]
## reorder factor levels of settings to match presentation in sec. 3.3 of the
# manuscript
experiments          <- arrange(experiments, setting)
experiments$setting  <- as.factor(experiments$setting)
experiments$form.int <- factor(experiments$form,
  labels=letters[seq_along(unique(experiments$form))])
experiments          <- transform(experiments,
  subsetting = interaction(setting, form.int))
# change levels of setting variable
cl <- lapply(cl, function(z) {
  z$setting <- z$subsetting <- z$form.int <- NULL
  z <- left_join(z, experiments)
})

## aggregated statistics for different settings (over all 500 simulation runs)
stats.df    <- do.call(rbind, lapply(cl, get_stats2))

rmean <- function(x) {
  round(mean(x), 3)
}
## LaTeX table with aggregated simulation results (presented in Appendix)
tab2 <- tabular(
  RowFactor(setting, name="Setting", spacing=1)*
  RowFactor(form.int, name="Model", spacing=1)*
  (RMSE+Coverage)*Heading()*rmean ~
  Factor(comparison, name="Comparison",
    levelnames=sub("Comparison ", "", levels(stats.df$comparison))),
  data=stats.df)
# output to tex file
tmp <- latex(tab2, file="simulation/simResultsv2.tex")


#### comparison plots and aggregated plots
## combine all setting in one data frame and add/improve labels
# for later use in graphics
df.cl <- do.call(rbind, cl)
df.cl <- transform(df.cl,
  below.ci = true.fit < fit - qnorm(0.975)*se,
  above.ci  = true.fit > fit + qnorm(0.975)*se)
df.cl$within.ci              <- !(df.cl$below.ci | df.cl$above.ci)
df.cl$ci.fit                 <- "within"
df.cl$ci.fit[df.cl$below.ci] <- "below"
df.cl$ci.fit[df.cl$above.ci] <- "above"
df.cl$ci.fit                 <- as.factor(df.cl$ci.fi)
df.cl$group                  <- "Data from Setting IV"
df.cl$group[df.cl$subsetting %in%
c("Setting I.a", "Setting I.b", "Setting II.a", "Setting II.b",
  "Setting III.c", "Setting III.d")] <- "Correct specification"
df.cl$group[df.cl$subsetting %in%
c("Setting I.c", "Setting I.d", "Setting II.c", "Setting II.d")] <- "Lead too long"
df.cl$group[df.cl$subsetting %in% c("Setting III.a", "Setting III.b")] <- "Lead too short"
df.cl$Penalty <- "P1"
df.cl$Penalty[df.cl$form.int %in% c("b", "d")] <- "P2"
df.cl$Lag <- "dynamic"
df.cl$Lag[df.cl$form.int %in% c("c", "d")] <- "static"
df.cl$Specification <- ifelse(df.cl$subsetting %in%
  paste("Setting", c("I.a", "I.b", "II.a", "II.b", "III.c", "III.d")),
  "\n correct", "\n incorrect")


## data frame with aggregated statistics for all settings (over all 500 simulation runs)
df.agg   <- get_stats(df.cl)
group.df <- unique(df.cl[, c("setting", "subsetting", "group", "Penalty", "Lag", "Specification")])
df.agg   <- left_join(df.agg, group.df)

#### comparison plots (containing individual trajectories)

## Data from Setting I
## correctly specified lag/lead, Penalties P1 and P2
Ia <- filter(df.cl, subsetting == "Setting I.a")
Ib <- filter(df.cl, subsetting == "Setting I.b")
gg.Ia <- ggcomparison_lines2(Ia, alpha=0.3)+
  stat_summary(fun.y="mean", geom="line", col="grey90", lty=2, lwd=1.2) +
  ggtitle(as.character(unique(Ia$subsetting)))
gg.Ib <- ggcomparison_lines2(Ib, alpha=0.3)+
  stat_summary(fun.y="mean", geom="line", col="grey90", lty=2, lwd=1.2) +
  ggtitle(as.character(unique(Ib$subsetting)))


jpeg("simulation/SettingIaIb.jpeg", width=600, height=800)
grid.draw(
  rbind(
    ggplotGrob(gg.Ia),
    ggplotGrob(gg.Ib),
    size = "last"))
dev.off()

pdf("simulation/SettingIaIb.pdf", width=8, height=13)
grid.draw(
  rbind(
    ggplotGrob(gg.Ia),
    ggplotGrob(gg.Ib),
    size = "last"))
dev.off()

## lead time too long, Penalties P1 and P2
Ic <- filter(df.cl, subsetting == "Setting I.c")
Id <- filter(df.cl, subsetting == "Setting I.d")
gg.Ic <- ggcomparison_lines2(Ic, alpha=0.3)+
  stat_summary(fun.y="mean", geom="line", col="grey90", lty=2, lwd=1.2)+
  ggtitle(as.character(unique(Ic$subsetting)))
gg.Id <- ggcomparison_lines2(Id, alpha=0.3)+
  stat_summary(fun.y="mean", geom="line", col="grey90", lty=2, lwd=1.2)+
  ggtitle(as.character(unique(Id$subsetting)))

pdf("simulation/SettingIcId.pdf", width=8, height=12)
grid.draw(
  rbind(
    ggplotGrob(gg.Ic),
    ggplotGrob(gg.Id),
    size = "last"))
dev.off()

jpeg("simulation/SettingIcId.jpeg", width=600, height=800)
grid.draw(
  rbind(
    ggplotGrob(gg.Ic),
    ggplotGrob(gg.Id),
    size = "last"))
dev.off()

## data from Setting II
## correctly specified lag/lead, Penalties P1 and P2
IIa <- filter(df.cl, subsetting == "Setting II.a")
IIb <- filter(df.cl, subsetting == "Setting II.b")
gg.IIa <- ggcomparison_lines2(IIa, alpha=0.3) +
  stat_summary(fun.y="mean", geom="line", col="grey90", lty=2, lwd=1.2)+
  ggtitle(as.character(unique(IIa$subsetting)))
gg.IIb <- ggcomparison_lines2(IIb, alpha=0.3) +
  stat_summary(fun.y="mean", geom="line", col="grey90", lty=2, lwd=1.2)+
  ggtitle(as.character(unique(IIb$subsetting)))

pdf("simulation/SettingIIaIIb.pdf", width=8, height=12)
grid.draw(
  rbind(
    ggplotGrob(gg.IIa),
    ggplotGrob(gg.IIb),
    size = "last"))
dev.off()

jpeg("simulation/SettingIIaIIb.jpeg", width=600, height=800)
grid.draw(
  rbind(
    ggplotGrob(gg.IIa),
    ggplotGrob(gg.IIb),
    size = "last"))
dev.off()

## lead too long, Penalties P1 and P2
IIc <- filter(df.cl, subsetting == "Setting II.c")
IId <- filter(df.cl, subsetting == "Setting II.d")
gg.IIc <- ggcomparison_lines2(IIc, alpha=0.3)+
  stat_summary(fun.y="mean", geom="line", col="grey90", lty=2, lwd=1.2)+
  ggtitle(as.character(unique(IIc$subsetting)))
gg.IId <- ggcomparison_lines2(IId, alpha=0.3)+
  stat_summary(fun.y="mean", geom="line", col="grey90", lty=2, lwd=1.2)+
  ggtitle(as.character(unique(IId$subsetting)))
pdf("simulation/SettingIIcIId.pdf", width=8, height=12)
grid.draw(
  rbind(
    ggplotGrob(gg.IIc),
    ggplotGrob(gg.IId),
    size = "last"))
dev.off()
jpeg("simulation/SettingIIcIId.jpeg", width=600, height=800)
grid.draw(
  rbind(
    ggplotGrob(gg.IIc),
    ggplotGrob(gg.IId),
    size = "last"))
dev.off()

## data from Setting III
## correctly specified lag/lead, Penalties P1 and P2
IIIa <- filter(df.cl, subsetting == "Setting III.a")
IIIb <- filter(df.cl, subsetting == "Setting III.b")
gg.IIIa <- ggcomparison_lines2(IIIa, alpha=0.3) +
  stat_summary(fun.y="mean", geom="line", col="grey90", lty=2, lwd=1.2)+
  ggtitle(as.character(unique(IIIa$subsetting)))
gg.IIIb <- ggcomparison_lines2(IIIb, alpha=0.3) +
  stat_summary(fun.y="mean", geom="line", col="grey90", lty=2, lwd=1.2)+
  ggtitle(as.character(unique(IIIb$subsetting)))

pdf("simulation/SettingIIIaIIIb.pdf", width=8, height=12)
grid.draw(
  rbind(
    ggplotGrob(gg.IIIa),
    ggplotGrob(gg.IIIb),
    size = "last"))
dev.off()
jpeg("simulation/SettingIIIaIIIb.jpeg", width=600, height=800)
grid.draw(
  rbind(
    ggplotGrob(gg.IIIa),
    ggplotGrob(gg.IIIb),
    size = "last"))
dev.off()

## lead too short, Penalties P1 and P2)
IIIc <- filter(df.cl, subsetting == "Setting III.c")
IIId <- filter(df.cl, subsetting == "Setting III.d")
gg.IIIc <- ggcomparison_lines2(IIIc, alpha=0.3)+
  stat_summary(fun.y="mean", geom="line", col="grey90", lty=2, lwd=1.2)+
  ggtitle(as.character(unique(IIIc$subsetting)))
gg.IIId <- ggcomparison_lines2(IIId, alpha=0.3)+
  stat_summary(fun.y="mean", geom="line", col="grey90", lty=2, lwd=1.2)+
  ggtitle(as.character(unique(IIId$subsetting)))
pdf("simulation/SettingIIIcIIId.pdf", width=8, height=12)
grid.draw(
  rbind(
    ggplotGrob(gg.IIIc),
    ggplotGrob(gg.IIId),
    size = "last"))
dev.off()

jpeg("simulation/SettingIIIcIIId.jpeg", width=600, height=800)
grid.draw(
  rbind(
    ggplotGrob(gg.IIIc),
    ggplotGrob(gg.IIId),
    size = "last"))
dev.off()

## data from null case (ELRA = 0)
# dynamic lag-lead, Penalties P1 and P2
IVa <- filter(df.cl, subsetting == "Setting IV.a")
IVb <- filter(df.cl, subsetting == "Setting IV.b")
gg.IVa <- ggcomparison_lines2(IVa, alpha=0.3) +
  stat_summary(fun.y="mean", geom="line", col="grey90", lty=2, lwd=1.2)+
  ggtitle(as.character(unique(IVa$subsetting)))
gg.IVb <- ggcomparison_lines2(IVb, alpha=0.3) +
  stat_summary(fun.y="mean", geom="line", col="grey90", lty=2, lwd=1.2)+
  ggtitle(as.character(unique(IVb$subsetting)))

pdf("simulation/SettingIVaIVb.pdf", width=8, height=12)
grid.draw(
  rbind(
    ggplotGrob(gg.IVa),
    ggplotGrob(gg.IVb),
    size = "last"))
dev.off()

jpeg("simulation/SettingIVaIVb.jpeg", width=600, height=800)
grid.draw(
  rbind(
    ggplotGrob(gg.IVa),
    ggplotGrob(gg.IVb),
    size = "last"))
dev.off()

# static lag-lead, Penalties P1 and P2
IVc <- filter(df.cl, subsetting == "Setting IV.c")
IVd <- filter(df.cl, subsetting == "Setting IV.d")
gg.IVc <- ggcomparison_lines2(IVc, alpha=0.3)+
  stat_summary(fun.y="mean", geom="line", col="grey90", lty=2, lwd=1.2)+
  ggtitle(as.character(unique(IVc$subsetting)))
gg.IVd <- ggcomparison_lines2(IVd, alpha=0.3)+
  stat_summary(fun.y="mean", geom="line", col="grey90", lty=2, lwd=1.2)+
  ggtitle(as.character(unique(IVd$subsetting)))

pdf("simulation/SettingIVcIVd.pdf", width=8, height=12)
grid.draw(
  rbind(
    ggplotGrob(gg.IVc),
    ggplotGrob(gg.IVd),
    size = "last"))
dev.off()

jpeg("simulation/SettingIVcIVd.jpeg", width=600, height=800)
grid.draw(
  rbind(
    ggplotGrob(gg.IVc),
    ggplotGrob(gg.IVd),
    size = "last"))
dev.off()

#### Overall summary graphics
df.agg$Comparison <- factor(df.agg$comparison, labels=LETTERS[1:6])

## Summarized RMSE results for each (Sub-)Setting
gg.rmse.agg <- ggplot(df.agg, aes(x=Comparison, y=RMSE)) +
    # geom_point(size=0.75) +
geom_line(aes(group=form.int, col=Penalty, lty=Lag), lwd=1.2) +
  scale_color_grey() +
  facet_wrap(~setting) +
  ggtitle("RMSE")
## Summarized Coverage results for each (Sub-)Setting
gg.cov.agg <- ggplot(df.agg, aes(x=Comparison, y=Coverage)) +
geom_hline(yintercept=0.95, lty=1, lwd=0.1) +
    # geom_point(size=0.75) +
geom_line(aes(group=form.int, col=Penalty, lty=Lag), lwd=1.2) +
  scale_color_grey() +
  facet_wrap(~setting) +
  ggtitle("Coverage") +
  ylim(c(min(df.agg$Coverage), 1))


pdf("simulation/simResAll.pdf", width=6, height=8)
grid.draw(
  rbind(
    ggplotGrob(gg.cov.agg),
    ggplotGrob(gg.rmse.agg),
    size = "last"))
dev.off()

# jpeg("simulation/simResAll.jpg")
# grid.draw(
#     rbind(
#         ggplotGrob(gg.cov.agg),
#         ggplotGrob(gg.rmse.agg),
#         size = "last"))
# dev.off()


## Overall summary graphics (grouped by subsetting)
df.agg$Comparison <- factor(df.agg$comparison, labels=LETTERS[1:6])
df.agg$Subsetting <- df.agg$form.int

## Summarized RMSE results for each (Sub-)Setting
gg.rmse.agg <- ggplot(df.agg, aes(x=Comparison, y=RMSE)) +
    # geom_point(size=0.75) +
geom_line(aes(group=form.int, lty=Lag, col=Penalty), lwd=1.2) +
  scale_color_grey() +
  facet_wrap(~setting, labeller=labeller(Specification=label_both)) +
  ggtitle("RMSE")
## Summarized Coverage results for each (Sub-)Setting
gg.cov.agg <- ggplot(df.agg, aes(x=Comparison, y=Coverage)) +
  geom_hline(yintercept=0.95, lty=1, lwd=0.3) +
  # geom_point(size=0.75) +
  geom_line(aes(group=form.int, col=Penalty, lty=Lag), lwd=1.2) +
  scale_color_grey() +
  facet_wrap(~setting, labeller=labeller(Specification=label_both)) +
  ggtitle("Coverage") +
  ylim(c(min(df.agg$Coverage), 1))


pdf("simulation/simResAll.pdf", width=6, height=8)
grid.draw(
  rbind(
    ggplotGrob(gg.cov.agg),
    ggplotGrob(gg.rmse.agg),
    size = "last"))
dev.off()

# jpeg("simulation/simResAllBySubsetting.jpg", width=550, height=450)
# grid.draw(
#     rbind(
#         ggplotGrob(gg.cov.agg),
#         ggplotGrob(gg.rmse.agg),
#         size = "last"))
# dev.off()



############################# Simulation Part B ################################
library(batchtools)
reg_dlnm <- loadRegistry("../../simulation/comparison/dlnm-tv-surv-registry/")

#### DLNM for survival data estimated by PAM
id_pam_ped <- findExperiments(prob.name="dlnm_sim_ped", algo.name="pam_dlnm_ped")
res_pam_ped <- reduceResultsDataTable(ids=findDone(id_pam_ped[,1])) %>%
  as_tibble() %>%
  tidyr::unnest()

stats_pam_ped <- res_pam_ped %>%
  mutate(
    mse = (fit-truth)^2,
    coverage = (truth <= fit + qnorm(0.975)*se) &
    (truth >= fit - qnorm(0.975)*se)) %>%
  group_by(job.id) %>%
  dplyr::summarize(
    RMSE = sqrt(mean(mse)),
    coverage = sum(coverage)/n()) %>%
  dplyr::summarize(
    RMSE = mean(RMSE),
    coverage = mean(coverage))
saveRDS(stats_pam_ped, "simulation/simstats_pam_dlnm.Rds")

av_pam_ped <- res_pam_ped %>%
  group_by(x, lag2) %>%
  summarize_at(vars(fit, truth), mean ) %>%
  arrange(x, lag2) %>%
  tidyr::gather(type, value, fit:truth)%>%
  mutate(type = case_when(
    type == "fit" ~ "PAM",
    TRUE ~ "TRUTH"))


gg_av_pam_ped <- ggplot(av_pam_ped, aes(x=x, y=lag2, fill=value)) +
  geom_tile() +
  geom_contour(aes(z=value), col="grey70") +
  scale_fill_gradient2(
    name = expression(h(t-t[e], z(t[e]))),
    low="firebrick", high="steelblue") +
  facet_wrap(~type) +
  ylab(expression(t-t[e])) + xlab(expression(z(t[e]))) +
  theme(
    axis.title   = element_text(size = rel(1.3)),
    axis.text    = element_text(size = rel(1.2)),
    legend.text  = element_text(size = rel(1.2)),
    legend.title = element_text(size = rel(1.3)),
    strip.text   = element_text(size = rel(1.3)))

ggsave("simulation/dlnm_ped.jpg", gg_av_pam_ped, device="jpeg", width=7, height=4)
ggsave("simulation/dlnm_ped.pdf", gg_av_pam_ped, device="pdf", width=7, height=4)


########################(time-varying) "complex" DLNM ##########################
reg_dlnm <- loadRegistry("../../simulation/comparison/dlnm-tv-surv-registry/")
#### time-varying DLNM for survival data estimated by PAM
id_pam_dlnm_ped_tv <- findExperiments(
  prob.name="sim_dlnm_ped_tv",
  algo.name="pam_dlnm_ped_tv")
res_pam_dlnm_ped_tv <- reduceResultsDataTable(ids=findDone(id_pam_dlnm_ped_tv[,1])) %>%
  as_tibble() %>%
  tidyr::unnest()

av_tv <- res_pam_dlnm_ped_tv %>%
  group_by(x, lag2, time_df) %>%
  summarize_at(vars(fit, truth), mean ) %>%
  arrange(time_df)

av_tv_sub <- filter(av_tv, time_df %in% c(1, 20, 40)) %>%
  tidyr::gather(type, value, fit, truth) %>%
  mutate(type = case_when(
    type == "fit" ~ "PAM",
    TRUE ~ "TRUTH"))

gg_av_tv <- ggplot(av_tv_sub, aes(x=x, y=lag2)) +
  geom_tile(aes(fill=value)) +
  geom_contour(aes(z = value), col="grey70") +
  scale_fill_gradient2(
    name = expression(h(t-t[e],z(t[e]))),
    low="firebrick", high="steelblue") +
  facet_grid(time_df ~ type) +
  ylab(expression(t-t[e])) + xlab(expression(z(t[e]))) +
  theme(
    axis.title   = element_text(size = rel(1.3)),
    axis.text    = element_text(size = rel(1.2)),
    legend.text  = element_text(size = rel(1.2)),
    legend.title = element_text(size = rel(1.3)),
    strip.text   = element_text(size = rel(1.3)))

ggsave("simulation/dlnm_ped_tv.pdf", gg_av_tv, width=8, height=10)

#### create animation
anim_df <- av_tv %>%
  mutate(
    group = row_number(),
    time  = round(time_df),
    time2 = factor(
      paste0("t = ", time),
      levels =  paste0("t = ", time)),
    ease  = "cubic-in-out") %>%
  select(-time_df) %>%
  tidyr::gather(type, value, fit:truth) %>%
  arrange(time) %>%
    mutate(type = case_when(
      type == "fit" ~ "PAM",
      TRUE ~ "TRUTH"))

gg_3d <- ggplot(anim_df, aes(x=x, y=lag2, fill=value, z=value, frame=time2)) +
  geom_tile() +
  geom_contour(aes(group=time), col="grey70") +
  scale_fill_gradient2(
    name = expression(h(t-t[e],z(t[e]))),
    low="steelblue", high="firebrick") +
  facet_wrap(~type) +
  ylab(expression(t-t[e])) + xlab(expression(z(t[e]))) +
  theme(
    legend.position = "bottom",
    axis.title   = element_text(size = rel(1.3)),
    axis.text    = element_text(size = rel(1.2)),
    legend.text  = element_text(size = rel(1.2)),
    legend.title = element_text(size = rel(1.3)),
    strip.text   = element_text(size = rel(1.3)))

# library(gganimate)
# gganimate(gg_3d, "simulation/fit-vs-truth-ped.gif",
#   interval=0.3, ani.width=700, height=300)