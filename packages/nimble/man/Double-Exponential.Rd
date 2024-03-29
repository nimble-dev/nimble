% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distributions_implementations.R
\name{Double-Exponential}
\alias{Double-Exponential}
\alias{ddexp}
\alias{rdexp}
\alias{pdexp}
\alias{qdexp}
\title{The Double Exponential (Laplace) Distribution}
\usage{
ddexp(x, location = 0, scale = 1, rate = 1/scale, log = FALSE)

rdexp(n, location = 0, scale = 1, rate = 1/scale)

pdexp(
  q,
  location = 0,
  scale = 1,
  rate = 1/scale,
  lower.tail = TRUE,
  log.p = FALSE
)

qdexp(
  p,
  location = 0,
  scale = 1,
  rate = 1/scale,
  lower.tail = TRUE,
  log.p = FALSE
)
}
\arguments{
\item{x}{vector of values.}

\item{location}{vector of location values.}

\item{scale}{vector of scale values.}

\item{rate}{vector of inverse scale values.}

\item{log}{logical; if TRUE, probability density is returned on the log scale.}

\item{n}{number of observations.}

\item{q}{vector of quantiles.}

\item{lower.tail}{logical; if TRUE (default) probabilities are \eqn{P[X \le x]}; otherwise, \eqn{P[X > x]}.}

\item{log.p}{logical; if TRUE, probabilities p are given by user as log(p).}

\item{p}{vector of probabilities.}
}
\value{
\code{ddexp} gives the density, \code{pdexp} gives the distribution
function, \code{qdexp} gives the quantile function, and \code{rdexp}
generates random deviates.
}
\description{
Density, distribution function, quantile function and random generation
for the double exponential distribution,
allowing non-zero location, \code{mu},
and non-unit scale, \code{sigma}, or non-unit rate, \code{tau}
}
\details{
See Gelman et al., Appendix A or the BUGS manual for mathematical details.
}
\examples{
x <- rdexp(50, location = 2, scale = 1)
ddexp(x, 2, 1)
}
\references{
Gelman, A., Carlin, J.B., Stern, H.S., and Rubin, D.B. (2004) \emph{Bayesian Data Analysis}, 2nd ed. Chapman and Hall/CRC.
}
\seealso{
\link{Distributions} for other standard distributions
}
\author{
Christopher Paciorek
}
