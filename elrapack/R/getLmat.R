#' Helper functions for creation of Lag-Lead matrices
#'
#' @rdname LagLead
#' @keywords internal
#' @export
createLmatDyn <- function(
  lead = 4,
  lag.c = 4,
  lag.f = 2,
  n.days = 11,
  brks = c(0:11, seq(15, 55, by = 5), 61.1), time.shift = 4)  {

  time.seq <- seq_len(n.days)
    if ( is.null(lag.c) ) {
        lag.c <- max(brks)
        lag.f <- 0
    }
  lag.vec  <- lag.f * time.seq + lag.c + time.seq
  lead.vec <- time.seq + lead

  lag.mat  <- sapply(lag.vec, function(z) z > brks[-length(brks)])
  lead.mat <- sapply(lead.vec, function(z) z <= brks[-1] )

  L        <- (lead.mat & lag.mat) * 1
  Lsub     <- L[-c(seq_len(time.shift), length(brks)), ]

  Lsub

}

#' @rdname LagLead
#' @keywords internal
#' @export
t_lead <- function(te, te.add.lead=0, lead.const=4, lead.factor=2) {
  lead.const + (te + te.add.lead)*lead.factor
}

#' @rdname LagLead
#' @keywords internal
#' @export
lag_lead_df <- function(
  te           = 1:12,
  te.add.lag   = 0,
  t.lag        = 4,
  te.add.lead  = 0,
  lead.const   = 4,
  lead.factor  = 2,
  interval.seq = c(0L:12L, seq(15, 55, by = 5), 61),
  labels       = TRUE) {

  lead    = t_lead(te, te.add.lead = te.add.lead, lead.const = lead.const,
    lead.factor = lead.factor)
  w.begin = (te + te.add.lag) + t.lag
  w.end   = get_end(w.begin, lead, max.end=max(interval.seq))

  data.frame(
    te        = te,
    lag       = t.lag,
    lead      = lead,
    w.begin   = w.begin,
    w.end     = w.end,
    int.begin = get_interval(w.begin, interval.seq=interval.seq, labels = labels),
    int.end   = get_interval(w.end, interval.seq=interval.seq, labels = labels)
    )
}

#' @rdname LagLead
#' @keywords internal
#' @export
get_end <- function(start, lead, max.end) {

  end <- start + lead
  end <- ifelse(end > max.end, max.end, end)
}

#' @rdname LagLead
#' @keywords internal
#' @export
get_interval <- function(
  x,
  interval.seq = c(0L:12L, seq(15, 55, by = 5), 61),
  labels       = FALSE) {


  ind <- cut(x, interval.seq, labels=FALSE)
  # if(any(is.na(ind))) ind[is.na(ind)] <- max(ind, na.rm=TRUE)
  if(labels) {
    ind <- int_info2(min.int=0, brks=interval.seq)[ind, "interval"]
  }

  ind

}

#' @rdname LagLead
#' @keywords internal
#' @export
int_info2 <- function(
  brks    = c(0:12, seq(15, 55, by=5), 61),
  min.int = 4) {

  intlen <- diff(brks)
  tstart <- c(0, cumsum(intlen)[-length(intlen)])
  tend   <- tstart + intlen

  tdf <- data.frame(
    tstart = tstart,
    tend   = tend,
    intlen = intlen)
  tdf          <- dplyr::mutate(tdf, intmid = tstart + intlen/2)
  tdf$interval <- paste0("(", tdf$tstart, ",", tdf$tend, "]")
  tdf$interval <- factor(tdf$interval, levels=tdf$interval, labels=tdf$interval)

  ind.keep <- which(tstart >= min.int)
  subset(tdf, tstart >= min.int)

}


#' creates one instance of Lag/Lead mat
#' @param te Numeric/Integer vector specifying the times at which exposure occurred.
#' @param te.add.lag A numeric constant added to te before application of lag time
#' @param t.lag A numeric constant, specifying the time (from \code{te}) before
#' \code{te} can affect hazard.
#' @param te.add.lead A numeric constant, added to te before application of lead time.
#' @param lead.const A numeric constant, specifying the constant amount of time
#' after \code{te + t.lag}, in which \code{te} can still affect hazard.
#' @param lead.factor If the lead time is dynamic, this factor can be set different
#' to zero, such that \code{t.lead=lead.const + lead.factor*te}.
#' @param interval.seq The break points dividing the follow up into intervals.
#' @param t.min If some intervals are not of interest only intervals for t > t.min are
#' returned.
#' @import checkmate dplyr
#' @return A data frame with intervals as first column and \code{length(te)}
#' columns specifying the lag/lead for each \code{te}.
#' @keywords internal
#' @export
create_Lmat <- function(
  te           = 1:12,
  te.add.lag   = 0,
  t.lag        = 4,
  te.add.lead  = 0,
  lead.const   = 4,
  lead.factor  = 2,
  interval.seq = c(0:12, seq(15, 55, by = 5), 61.1),
  t.min        = 0) {

  assert_integer(te, lower=1, any.missing=FALSE, unique=TRUE)
  assert_numeric(interval.seq, lower=0, any.missing=FALSE, min.len=2)
  assert_number(te.add.lag, lower=0, upper=max(interval.seq), finite=TRUE)
  assert_number(t.lag, lower=1, upper=max(interval.seq), finite=TRUE)
  assert_number(lead.const, lower=0, upper=max(interval.seq), finite=TRUE)
  assert_number(lead.factor, lower=0, upper=max(interval.seq), finite=TRUE)
  assert_number(t.min, lower=0, upper=max(interval.seq), finite=TRUE)

  # create lag-lead information matrix
  ldf <- lag_lead_df(te=te, te.add.lag=te.add.lag, te.add.lead=te.add.lead,
    t.lag=t.lag, lead.const=lead.const, lead.factor=lead.factor,
    interval.seq=interval.seq)

  ind.begin <- get_interval(ldf$w.begin, interval.seq=interval.seq)
  ind.end   <- get_interval(ldf$w.end, interval.seq=interval.seq)

  int.info  <- int_info2(brks=interval.seq, min.int=0)
  int.keep  <- int.info$interval[which(int.info$tstart >= t.min)]

  ints <- apply(cbind(ind.begin, ind.end), 1, function(z) {
    z.i <- int.info$interval[z[1]:z[2]]
    int.info$interval %in% z.i
  }) * 1

  ints <- data.frame(intsL=int.info$interval, Lcols=ints)
  filter(ints, intsL %in% int.keep)


}


##' creates Lag/Lead matrix for all observations in a data set

#' @param data The complete data set.
#' @param Lmat A data frame where first column specifies intervals and the other
# columns the Lag/Lead structure for each day of exposure.
#' @param merge.col.data The name of the column (in data) on which data and Lmat
# should be merged. Defaults to \code{"int"}.
#' @param merge.col.Lmat The name of the column (in Lmat) on which data and Lmat
#' should be merged. Defaults to \code{"intsL"}
#' @details:
#' Lmat, that only contains information for all unique intervals and all days of
#' exposure, is merged with data, such that correct rows of Lmat are added to the
#' data, but only the Lmat matrix (now with as many rows as data) is returned.
#' @importFrom stats setNames
#' @importFrom dplyr left_join
#' @import checkmate
#' @keywords internal
#' @export
Lmat_data <- function(
  data,
  Lmat,
  merge.col.data = "int",
  merge.col.Lmat = "intsL") {

  # check inputs
  assert_data_frame(data)
  assert_data_frame(Lmat)
  assert_set_equal(levels(data[[merge.col.data]]), levels(Lmat[[merge.col.Lmat]]))

  nrow.data <- nrow(data)
  int <- data["int"]
  rm(data)
  int <- left_join(int, Lmat, by=setNames(merge.col.Lmat, merge.col.data))

  if(!(nrow(int)==nrow.data)) {
    stop("Left join not successful, number of rows produced no equal to number of
      rows in data")
  }

  # return int, -1 because we only want the Lmat, not the intervals
  # this is necessary because we need to store L as a (numeric) matrix in the
  # data for it to be processed properly by mgcv::gam
  as.matrix(int[,-1])

}


#'  heatmap of the effect terms for relevant L-Areas
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#' @keywords internal
#' @export
heatAdequacy <- function(
  hm         = NULL,
  day.names  = 1:11,
  int.names  = NULL,
  high.col   = "steelblue",
  low.col    = "white",
  grid.col   = "lightgrey",
  title.char = "") {

  ## hm: heat matrix containing the values to be ploted via heatmap()
  m.hm <- melt(hm)
  m.hm$Var1 <- factor(m.hm$Var1, labels = int.names)
  m.hm$Var2 <- factor(m.hm$Var2, labels = day.names)

  ggplot(m.hm, aes(x = Var2, y = rev(Var1))) +
    geom_tile(aes(fill = value), colour = grid.col) +
    scale_fill_gradient(low = low.col, high = high.col) +
    xlab("Protocol day") +
    scale_y_discrete("Interval j", labels = rev(int.names)) +
    theme(legend.position = "none") +
    labs(title = title.char)

}
