## packages
library(dlnm)
library(tsModel)
library(dplyr)
library(ggplot2)

#### Functions for simulations
## mostly from
## http://www.ag-myresearch.com/2017_gasparrini_biomet.html

## helper functions for functional shape
fdnorm <- function(x) (dnorm(x,1.5,2)+1.5*dnorm(x,7.5,1))
wpeak2 <- function(lag) 15*dnorm(lag,8,10)
wdnorm <- function(lag) 5*(dnorm(lag,4,6)+dnorm(lag,25,4))
fcomplex <- function(x,lag) 0.8*(fdnorm(x)-fdnorm(5)) *
    ifelse(is.na(x),NA,ifelse(x>=5,wdnorm(lag), wpeak2(lag)))

## additional additive, non-linear component of the log-baseline hazard
f0 <- function(t) {
  t <- (t) / max(t) * 3/2*pi
  (-1*sin(t))-0.5
}

# FUNCTION TO COMPUTE THE CUMULATIVE EFFECT GIVEN AN EXPOSURE HISTORY
fcumeff <- function(hist,lag,fun) sum(do.call(fun,list(hist,lag)))

## save static parts of simulation in input folder
dir.create("input/", showWarnings=FALSE)

#### static dlnm-ped
## reuse simulation from
## http://www.ag-myresearch.com/2017_gasparrini_biomet.html

# REAL TEMPERATURE SERIES, STANDARDIZED IN 0-10
# chicagoNMMAPS data from dlnm package
xz <- chicagoNMMAPS$temp
xz <- (xz-min(xz))/diff(range(xz))*10

## set parameters
# number of exposures per subject
nz <- 80 # t_e = (1:80 -41)

# number of subjects
n <- 2000 # i = 1:n

# sample exposures from the empirical distribution of x
x <- sample(xz, nz*n, replace=T)

## id and time variables
time <- rep(1:nz,n)
id   <- rep(1:n, each=nz)
## calculate cumulative effect for all time-points and all subjects (= eta)
# each element of x_list contains exposure history with nz elements
x_list   <- split(x, id)
# each element of lag_list contains a Lag matrix, with lags 0 to 40
# lag = 1 means exposure was recorded one day before time at which we want to
# model the hazard, i.e. t-te = 1
# lag = 40 means exposure was recorded 40 days before the time at which we want
# to model the hazard, i.e. t-te = 40
lag_list <- lapply(x_list, Lag, 0:40)

# For each subject, go through each row of the Lag matrix and calculate
# the cumulative effect for each lag.
# Note: this will be NA for the fist 39 rows, as for the DLNM we need
# max(lag) observations at each time-point at which we want to model the hazard.
# The cumulative effec is the sum over all 40 partial effects, i.e.
# h(t-te=1, z[te=t-1]) + h(t-te=2, z[te=t-2]) + ...+ h(t-te=40, z[te=t-40])
eff_vec  <- sapply(lag_list, apply, 1, fcumeff, 0:40, "fcomplex")

## create full data set (observations for each subject and time-point)
Xdf <- data.frame(id = id, time=time, eta_dlnm=as.vector(eff_vec))
## create covariates for linear functionals
# these need to be matrices with repeated entries, for mgcv to recognize the
# specified terms as cumulative effects
# - Z = exposure history matrix
# - te_df = matrix of exposure times t_e
# - time_df = matrix of time-points of the follow-up
# - LL = lag-lead matrix
Xdf$Z       <- do.call(rbind, lapply(x_list, matrix, nrow=nz, ncol=nz, byrow=TRUE))
Xdf$te_df   <- matrix(1:nz, nrow=nrow(Xdf), ncol=nz, byrow=TRUE)
Xdf$time_df <- matrix(Xdf$time, nrow=nrow(Xdf), ncol=nz)
diff_df     <- (Xdf$time_df - Xdf$te_df)
Xdf$LL      <- ((diff_df >=0) & (diff_df <=40))*1
## create crossbasis for DLNM estimation by Gasparrini
lag_mat <- do.call(rbind, lag_list)
cb <- crossbasis(lag_mat[,1], lag=c(0,40),
	argvar = list(fun='ps',df=9), arglag=list(fun='ps'))
Xdf$cb      <- cb

## follow up starts after 40 days of exposure, such that every subject has
# complete exposure history of 40 exposures at the beginning of the follow-up
Xdf$time <- Xdf$time - 41
## remove all obs before follow-up starts
Xdf        <- Xdf[Xdf$time >= 0, ]
Xdf$tstart <- Xdf$time
Xdf$tend   <- Xdf$tstart + 1
# calculate linear predictor using cumulative effect + intercept + smooth baseline
Xdf$eta  <- -3.5 + f0(Xdf$time) + Xdf$eta_dlnm

#### create new data for prediction + truth column
ndf <- expand.grid(Z = seq(0,10,0.25), te_df = 0:40)
# any time will do (later we only use coefficients for f_complex estimation)
ndf$tend  <- 20
ndf$LL    <- 1
ndf$truth <- apply(ndf, 1, function(x) {fcomplex(x[1], x[2])})
# ggplot(ndf, aes(x=Z, y=te_df)) +
# 	geom_tile(aes(fill=truth)) +
# 	geom_contour(aes(z=truth)) +
# 	scale_fill_gradient2(low="firebrick", high="steelblue")

#### dlnm stuff
cb_attr <- attributes(cb)

dlnm_ped_static <- list(
	X        = Xdf,
	ndf      = ndf,
	cb_attr  = cb_attr,
	cbgamPem = cbPen(cb))

saveRDS(dlnm_ped_static, "input/static_dlnm_ped.Rds")



#### static for time-varying dlnm for survival data
## helper functions for functional shape
fcumeff_tv <- function(hist,lag,fun, t, tmax)
	sum(do.call(fun,list(hist, lag, t, tmax)))

## functions for time-varying dlnm
ft <- function(t, tmax) {
	-1*cos(t/tmax*pi)
}
## 3d, time-varying complex surface function
fcomplex_tv <- function(x,lag, t, tmax) {
	0.8*(fdnorm(x)-fdnorm(5)) *
    ifelse(is.na(x),NA,ifelse(x>=5,wdnorm(lag), wpeak2(lag)))*
    ft(t, tmax)
  }

f0 <- function(t) {
  t <- (t) / max(t) * 3/2*pi
  (-.5*sin(t))
}
# REAL TEMPERATURE SERIES, STANDARDIZED IN 0-10
xz <- chicagoNMMAPS$temp
xz <- (xz-min(xz))/diff(range(xz))*10

## set parameters
# number of exposures
nz <- 80 # t_e = 1:80
# number of subjects
n <- 5000 # i = 1:n
# sample exposures from the empirical distribution of x
x <- sample(xz, nz*n, replace=T)
## id and time variables
time <- rep(1:nz,n)
id   <- rep(1:n, each=nz)
## calculate cumulative effect for all time-points and all subjects (= eta)
x_list   <- split(x, id)
lag_list <- lapply(x_list, Lag, 0:40)
eff_vec  <- sapply(lag_list,
	function(Q) {
		ceff <- sapply(seq_len(nrow(Q)), function(i) {
			fcumeff_tv(hist=Q[i,], lag=0:40, fun="fcomplex_tv", t=i-40, tmax=40)
	})
		ceff
	})

## create full data set (observations for each subject and time-point)
Xdf <- data.frame(id = id, time = time, eta_dlnm = as.vector(eff_vec))
## create covariates for linear functionals
# - Z = exposure history matrix
# - te_df = matrix of exposure time t_e
# - time_df = matrix of time-points of the follow-up
# - LL = lag-lead matrix
Xdf$Z       <- do.call(rbind, lapply(x_list, matrix, nrow=nz, ncol=nz, byrow=TRUE))
Xdf$te_df   <- matrix(1:nz, nrow=nrow(Xdf), ncol=nz, byrow=TRUE)
Xdf$time_df <- matrix(Xdf$time, nrow=nrow(Xdf), ncol=nz)
diff_df     <- (Xdf$time_df - Xdf$te_df)
Xdf$LL      <- ((diff_df >=0) & (diff_df <=40))*1
## follow up starts after 40 days of exposure, such that every subject has
# complete exposure history of 40 exposures at the beginning of the follow-up
Xdf$time <- Xdf$time - 41
## remove all obs before follow-up starts
Xdf         <- Xdf[Xdf$time >= 0, ]
Xdf$tstart  <- Xdf$time
Xdf$tend    <- Xdf$tstart + 1
# calculafe linear predictor using cumulative effect + intercept
Xdf$eta     <- -3 + f0(Xdf$time) + Xdf$eta_dlnm
Xdf$te_df   <- (Xdf$time_df - Xdf$te_df)*Xdf$LL
Xdf$time_df <- Xdf$time_df - 40

#### create new data for prediction + truth column
ndf <- expand.grid(Z = seq(0,10,0.25), te_df = 0:40, time_df=1:40)
# any time will do (later we only use coefficients for f_complex estimation)
ndf$tend  <- 20
ndf$LL <-1
ndf$truth <- apply(ndf, 1, function(x) {fcomplex_tv(x[1], x[2], x[3], tmax=40)})

dlnm_ped_static_tv <- list(
	X   = Xdf,
	ndf = ndf)

saveRDS(dlnm_ped_static_tv, "input/static_dlnm_ped_tv.Rds")