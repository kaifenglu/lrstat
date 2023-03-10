% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wrappers.R
\name{qtpwexp}
\alias{qtpwexp}
\title{Quantile function of truncated piecewise exponential distribution}
\usage{
qtpwexp(
  probability,
  piecewiseSurvivalTime = 0,
  lambda = 0.0578,
  lowerBound = 0
)
}
\arguments{
\item{probability}{The scalar probability corresponding to the quantile.}

\item{piecewiseSurvivalTime}{A vector that specifies the time intervals for
the piecewise exponential survival distribution. Must start with 0, e.g.,
\code{c(0, 6)} breaks the time axis into 2 event intervals:
[0, 6) and [6, Inf).
Defaults to 0 for exponential distribution.}

\item{lambda}{A vector of hazard rates for the event. One for
each analysis time interval.}

\item{lowerBound}{The left truncation time point for the survival time.
Defaults to 0 for no truncation.}
}
\value{
The quantile x such that
P(X > x | X > lowerBound) = 1 - probability.
}
\description{
Obtains the quantile of a piecewise expoenential distribution
given that it exceeds a specified lower bound.
}
\examples{
qtpwexp(probability = c(0.3, 0.5), piecewiseSurvivalTime = c(0, 6, 9, 15),
        lambda = c(0.025, 0.04, 0.015, 0.007))

}
