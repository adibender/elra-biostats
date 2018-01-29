#' Plot Lag-Lead windows
#'
#' Given a matrix defining a lag lead window, returns respective plot as a
#' \code{ggplot2} object.
#'
#' @param hm A lag-lead matrix
#' @param day.names Labels applied to exposure times
#' @param int.names Labels applied to intervals in which follow-up was partitioned
#' @param high.col Color used to highlight exposure times within the lag-lead window
#' @param low.col Color of exposure times outside the lag-lead window.
#' @param grid.col Color of grid lines.
#' @param title.char Will be added at plot title
#' @import ggplot2
#' @importFrom reshape2 melt
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
