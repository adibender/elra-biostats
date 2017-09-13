 sim_dlnm_ped <- function(
  data,
  job) {

  X <- data$X

  sdf      <- data.frame(id = X$id, rate = exp(X$eta), t = X$time)
  rt.l     <- split(sdf[,-1], f=sdf$id)
  rt.l     <- lapply(rt.l, as.list)
  new_time <- vapply(rt.l, do.call, what=msm::rpexp, numeric(1))

  sv_df <- data.frame(id=seq_len(length(rt.l)), surv = new_time)
  X <- merge(X, sv_df)
  X            <- X[X$time >=0, ]
  X$offset     <- 0
  X$ped_status <- 0
  X                     <- X[X$tstart < X$surv,]
  last_id               <- cumsum(rle(X$id)$lengths)
  X$ped_status[last_id] <- 1*(X$surv[last_id] < 40)
  X$offset[last_id]     <- log(X$surv[last_id] - X$tstart[last_id])
  X$te_df               <- (X$time_df - X$te_df)*X$LL

  return(X)

}


## as sim_dlnm_wrapper, but for survival times
 sim_dlnm_ped_tv <- function(
  data,
  job) {

  X <- data$X

  sdf      <- data.frame(id = X$id, rate = exp(X$eta), t = X$time)
  rt.l     <- split(sdf[,-1], f=sdf$id)
  rt.l     <- lapply(rt.l, as.list)
  new_time <- vapply(rt.l, do.call, what=msm::rpexp, numeric(1))

  sv_df <- data.frame(id=seq_len(length(rt.l)), surv = new_time)
  X <- merge(X, sv_df)
  X            <- X[X$time >=0, ]
  X$offset     <- 0
  X$ped_status <- 0
  X                     <- X[X$tstart < X$surv,]
  last_id               <- cumsum(rle(X$id)$lengths)
  X$ped_status[last_id] <- 1*(X$surv[last_id] < max(X$time_df))
  X$offset[last_id]     <- log(X$surv[last_id] - X$tstart[last_id])

  return(X)

}