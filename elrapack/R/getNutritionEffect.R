#' Extract effect differences for 2 hypothetical patients
#'
#' @param effectname Name of variable(s) for which the effect should be extracted.
#' @import mgcv
#' @importFrom stats predict
#' @keywords internal
#' @export
getEffectDiffs <- function(
  m,
  newdata1,
  newdata2,
  effectname = "AdequacyCals") {

  ## deal with fits created with old mgcv before ti() was implemented
  ## when loaded with newer versions that have ti(): add the missing "mc"-slot
  if(!is.null(mgcv::ti)){
    m$smooth <- lapply(m$smooth, function(trm){
      if("margin" %in% names(trm) & is.null(trm$mc)){
        trm$mc <- rep(TRUE, length(trm$margin))
                #message("added mc slot to ", trm$label)
      }
      trm
    } )
  }

  useColumns <- grep(effectname, names(m$coefficients))
  covCoefs <- m$Vp[useColumns, useColumns]
  X1 <- predict(m, newdata = newdata1, type = "lpmatrix")[,useColumns]
  X2 <- predict(m, newdata = newdata2, type = "lpmatrix")[,useColumns]
  X <- X2 - X1
  fit <- drop(X %*% m$coefficients[useColumns])
  se <- sqrt(rowSums((X %*% covCoefs) * X))
  hi <- fit + 2 * se
  lo <- fit - 2 * se
  cbind(intmid = newdata1$intmid, fit = fit, se = se, hi = hi, lo = lo)

}
