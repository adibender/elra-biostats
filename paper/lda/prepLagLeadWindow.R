library(grid)
library(gridExtra)
library(ggplot2)
theme_set(theme_bw())

nutri        <- readRDS("../../dataCurrent/nutriOrigSmall.Rds")
seq.ints     <- seq_along(unique(nutri$intmid))
colindex     <- if(ncol(nutri$LHartlDynf)==12) -1 else 1:11
## graphic titles
dyn.title <- expression(
  atop("Dynamic time window",
    atop("lag = 4 days; lead = 4 + 2 x n days",
      "(n = number of days in the ICU)")))
static.title <- expression(atop("Static time window",
  atop("lag = 4 days; lead = 30 days")))

pdf("lagLeadWindow.pdf", width=7, height=6)
grid.arrange(
  heatAdequacy(
    nutri$LHartlDynf[seq.ints, colindex],
    high.col   = "grey20",
    title.char = dyn.title, day.names = 1:11,
    int.names  = unique(elrapack::int_info(brks = c(0:30))$interval)) +
   xlab(expression(paste("Nutrition day ", t[e]))) +
   theme(
      axis.text  = element_text(size = rel(1.2)),
      axis.title = element_text(size = rel(1.3))),
  heatAdequacy(nutri$LSimple460f[seq.ints, colindex],
    high.col   = "grey20",
    title.char = static.title, day.names = 1:11,
    int.names  = unique(elrapack::int_info(brks = c(0:30))$interval)) +
  xlab(expression(paste("Nutrition day ", t[e]))) +
  theme(
      axis.text  = element_text(size = rel(1.2)),
      axis.title = element_text(size = rel(1.3))),
  nrow = 1)
dev.off()