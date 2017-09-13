library(mgcv)
library(ggplot2)
library(grid)
library(gridExtra)
theme_set(theme_bw())

static <- readRDS("../../simulation/modelEvaluation/input/static.Rds")

coef.SettingI <- static$coefs.list[[1]]

X         <- static$Xcompare.list[[1]][[1]]
f0.df     <- data.frame(intmid=4.5:29.5)
ind.f0    <- grep("s(intmid)", names(coef.SettingI), fixed=TRUE)
f0.df$eff <- X[,ind.f0] %*% coef.SettingI[ind.f0]

gg.f0 <- ggplot(f0.df, aes(x=intmid, y=eff)) +
	geom_line() +
	xlab("t") + ylab(expression(log(lambda[0](t))))

x.range <- range(static$complete.data$ApacheIIScore)
x.vec   <- seq(x.range[1], x.range[2], by=.5)

fx.df <- data.frame(
	intmid = rep(f0.df$intmid, times = length(x.vec)),
	x      = rep(x.vec, each = nrow(f0.df)))
## could be anything, just needed to create model matrix
fx.df$y <- rnorm(nrow(fx.df))

mm.fx <- model.matrix(gam(y ~ s(intmid, bs="ps") + te(x, intmid), data=fx.df))

## columns of design matrix and indices of vector for elements of f(x,t)
ind.mm    <- grep("te(x,intmid)", colnames(mm.fx), fixed=TRUE)
ind.x     <- grep("te(ApacheIIScore,intmid)", names(coef.SettingI), fixed=TRUE)
fx.df$eff <- mm.fx[, ind.mm] %*% coef.SettingI[ind.x]

gg.fx <- ggplot(fx.df, aes(x=intmid, y=x, z=eff)) +
	geom_raster(aes(fill=eff)) +
	scale_fill_gradient2(
		name=expression(f(list(x,t))),
		low="steelblue", high="firebrick2") +
	geom_contour(col="grey30") +
	xlab("t") +
	ylab("x")


pdf("simulation/trueSimBaseConfounder.pdf", width=8, height=4)
grid.draw(
    cbind(
        ggplotGrob(gg.f0),
        ggplotGrob(gg.fx),
        size = "first"))
dev.off()


#### true ELRA effects
td.expertt1 <- get_true_diffs(static, "Expertt1")
td.expertt2 <- get_true_diffs(static, "Expertt2")
td.statict1 <- get_true_diffs(static, "Statict1")


gg.e1 <- ggplot(td.expertt1, aes(x=intmid, y=true.fit)) +
	geom_line() +
	facet_wrap(~label) +
	ylim(c(-2, 2)) +
	xlab("t") +
	ylab(expression(e[j])) +
	geom_hline(yintercept=0, lty=2) + ggtitle("Setting I") +
  scale_y_continuous(
    breaks = log(sy <- c(0.25, 0.5, 0.75, 1, 1.25, 2, 4)),
    labels = sy,
    limits = c(-1.5, 1.5))


gg.e2 <- ggplot(td.expertt2, aes(x=intmid, y=true.fit)) +
	geom_line() +
	facet_wrap(~label) +
	ylim(c(-2, 2)) +
	xlab("t") +
	ylab(expression(e[j])) +
	geom_hline(yintercept=0, lty=2) + ggtitle("Setting II") +
  scale_y_continuous(
    breaks = log(sy <- c(0.25, 0.5, 0.75, 1, 1.25, 2, 4)),
    labels = sy,
    limits = c(-1.5, 1.5))


gg.s1 <- ggplot(td.statict1, aes(x=intmid, y=true.fit)) +
	geom_line() +
	facet_wrap(~label) +
	ylim(c(-2, 2)) +
	xlab("t") +
	ylab(expression(e[j])) +
	geom_hline(yintercept=0, lty=2) + ggtitle("Setting III") +
  scale_y_continuous(
    breaks = log(sy <- c(0.25, 0.5, 0.75, 1, 1.25, 2, 4)),
    labels = sy,
    limits = c(-1.5, 1.5))


pdf("simulation/trueELRAall.pdf", width=8, height=13)
grid.draw(
    rbind(
        ggplotGrob(gg.e1),
        ggplotGrob(gg.e2),
        ggplotGrob(gg.s1),
        size = "last"))
dev.off()
