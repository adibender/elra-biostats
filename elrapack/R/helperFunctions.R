#' create start/end times, given breaks
#'
#' @import checkmate
#' @importFrom dplyr mutate filter
#' @keywords internal
#' @export
int_info <- function(
  brks    = c(0:12, seq(15, 55, by=5), 61),
  min.int = 4L) {

  # check inputs
  assert_numeric(brks, lower = 0, any.missing = FALSE)
  assert_int(min.int, lower  = 0L)

  intlen <- diff(brks)
  tstart <- brks[-length(brks)]
  tend   <- tstart + intlen

  tdf <- data.frame(
    tstart = tstart,
    tend   = tend,
    intlen = intlen)
  tdf <- mutate(tdf, intmid = tstart + intlen/2)
  tdf$interval <- paste0("(", tdf$tstart, ",", tdf$tend, "]")

  filter(tdf, tstart >= min.int)

}

#' get index of last observation of duplicated values
#'
#' @importFrom checkmate assert_integer
#' @keywords internal
#' @export
get_last <- function(id) {

  assert_integer(id, lower=0, any.missing=FALSE, all.missing=FALSE, min.len=2)
  if(is.unsorted(id)) stop("Input variable must be sorted!")

  cumsum(rle(id)$lengths)

}

