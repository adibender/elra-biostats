#### DLNM PED
dlnm_dlnm_ped <- function(
	job,
	data,
	instance,
	qn = qnorm(0.975),
	debug = FALSE) {

	require(mgcv)
	require(dlnm)
	require(dplyr)

	cbgamPen <- data$cbgamPen
	instance <- na.omit(instance)
	# form <- as.formula(form)
	mod      <- gam(ped_status ~ s(tend, k=6) + cb, data=instance,
		family = poisson(), offset=offset, method="REML", paraPen=list(cb=cbgamPen),
		control = list(trace=TRUE))

	cb <- instance$cb
	attr_ind <- !(names(data$cb_attr) %in% c("dim", "dimnames"))
	for(i in names(data$cb_attr[attr_ind])) {
		attr(cb, i) <- data$cb_attr[[i]]
	}

	cp       <- crosspred(cb,mod,from=0,to=10,by=0.25,bylag=1,cen=5)

	# STORE THE RESULTS
	ndf <- data$ndf[c(1,2,5)]
	colnames(ndf) <- c("x", "lag", "truth")
	ndf$fit <- reshape2::melt(cp$matfit)$value
	ndf$se <- reshape2::melt(cp$matse)$value
	ndf <- mutate(ndf,
		bias = fit - truth,
		cov  = (truth >= fit -qn*se & truth <= fit +qn*se),
		rmse = (fit - truth)^2)

	if(debug) {
		return(list(mod=mod, res=ndf))
	} else {
		return(ndf)
	}

}



pam_dlnm_ped <- function(
	job,
	data,
	instance,
	cen   = 5,
	debug = FALSE) {

	require(mgcv)
	require(dplyr)

	mod <- gam(
		ped_status ~ s(tend, k=6) + ti(Z, te_df, by=LL, k=c(9, 9), mc=c(1,0)),
		data=instance, family=poisson(), offset=offset, control = list(trace=TRUE))

	ndf    <- ndf2 <- data$ndf
	ndf2$Z <- cen
	X1     <- predict(mod, newdata=ndf, type="lpmatrix")[, -1]
	X      <- X1
	X2     <- predict(mod, newdata=ndf2, type="lpmatrix")[, -1]
	X      <- X1 - X2
	ndf    <- ndf %>% select(Z, te_df, LL, everything())
	colnames(ndf)[1:3] <- c("x", "lag", "LL")
	ndf$fit <- as.numeric(X %*% coef(mod)[-1])
	ndf$se  <- sqrt(rowSums((X %*% mod$Vp[-1,-1]) * X))
	ndf     <- ndf %>%  mutate(
		lag2 = lag,
		lag  = paste0("lag", lag))

	gc()

	if(debug) {
		return(list(mod = mod, ndf=ndf))
	} else {
		return(ndf)
	}

}

pam_dlnm_ped_tv <- function(
	job,
	data,
	instance,
	cen   = 5,
	debug = FALSE) {

	require(mgcv)
	require(dplyr)
	system.time({
		mod <- gam(
			ped_status ~s(tend, k=6) +
			ti(Z, te_df, time_df, by=LL, k=c(9,9,4), mc=c(1,0, 1)),
			data=instance, family=poisson(), offset=offset,
			control = list(trace=TRUE), method="REML")
	})
	ndf      <- ndf2      <- data$ndf
	ndf$tend <- ndf2$tend <- ndf$time_df
	ndf2$Z   <- cen
	X1  <- predict(mod, newdata=ndf, type="lpmatrix")[, -1]
	X   <- X1
	X2  <- predict(mod, newdata=ndf2, type="lpmatrix")[, -1]
	X   <- X1 - X2
	ndf <- ndf %>% select(Z, te_df, LL, everything())
	colnames(ndf)[1:3] <- c("x", "lag", "LL")
	ndf$fit <- as.numeric(X %*% coef(mod)[-1])
	ndf$se  <- sqrt(rowSums((X %*% mod$Vp[-1,-1]) * X))
	ndf <- ndf %>%  mutate(
		lag2 = lag,
		lag  = paste0("lag", lag))

	gc()

	if(debug) {
		return(list(mod = mod, ndf=ndf))
	} else {
		return(ndf)
	}

}