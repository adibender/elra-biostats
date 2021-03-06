% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eval_mgcv_teststat.R
\name{testEffectDiffs}
\alias{testEffectDiffs}
\title{Perform p-value calculation for difference between two hypothetical patients.}
\usage{
testEffectDiffs(m, newdata1, newdata2, effectname = "AdequacyCals",
  debug = FALSE, mod.df = FALSE)
}
\arguments{
\item{effectname}{Name of variable(s) for which the effect should be extracted.
\code{"AdequacyCals"}.}
}
\description{
Perform p-value calculation for difference between two hypothetical patients.
}
\keyword{internal}
