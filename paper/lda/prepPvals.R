## 
library(ggplot2)
theme_set(theme_bw())
library(magrittr)
library(tidyr)
library(dplyr)


####  Evaluate p-value simulation 


## Hypothesis Testing Overall Test (built-in):
pvals.standard <- readRDS("simulation/pvalsStandard.Rds")
std.pval.df    <- do.call(rbind, pvals.standard)

std.pval.df %<>% 
	select(-algo, -truth, -.count, -prob, -subsetting) %>% 
	rename(Subsetting=form.int) %>% 
	gather(Effect, pval, CII, CIII)

std.pval.df$Penalty <- ifelse(std.pval.df$Subsetting %in% c("a", "c"), "P1", "P2")
std.pval.df$Lag     <- ifelse(std.pval.df$Subsetting %in% c("a", "b"), "dynamic", "static")

std.agg.df <- std.pval.df %>% 
	group_by(setting, Penalty, Lag, Effect) %>% 
	summarise(
		rejected = round(mean(pval < 0.05), 3)) %>%
	mutate(
		x = 0.01,
		y = ifelse(Effect=="CII", 0.9, 0.8))

pdf("simulation/pvalsStandard.pdf", width=7.2, height=7.2)
# jpeg("simulation/pvalsStandard.jpg", width=600, height=600)
ggplot(filter(std.pval.df)) + 
	geom_abline(slope=1, intercept=0, linetype=2) + 
	geom_qq(aes(sample=pval, pch=Effect, col=Effect), size=0.9,
		distribution=qunif) + 
	facet_grid(Penalty + Lag ~ setting, 
		labeller=labeller(Penalty=label_both, Lag=label_both))+
	xlab("Expected") + ylab("Observed") + 
	scale_shape_manual(values=c(3, 19)) + 
	scale_color_manual(values=c("black", "grey70")) + 
	geom_text(
		data=filter(std.agg.df), 
		aes(x=x, y=y,
			label=paste("list(", Effect, ": alpha[e] ==", rejected, ")")), 
		parse=TRUE, size=3, hjust=0) + 
	scale_fill_manual(values=c("black", "grey70")) + 
	theme(axis.text.x = element_text(angle=90, hjust=1)) + 
	geom_hline(yintercept=0.05, lty=2, col="grey70") + 
	guides(colour=guide_legend(override.aes=list(alpha=1))) + 
	theme(legend.position="bottom")
dev.off()