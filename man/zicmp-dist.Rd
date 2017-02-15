\name{ZICMP Distribution}
\alias{rzicmp}
\alias{dzicmp}
\alias{pzicmp}
\alias{qzicmp}

\title{ZICMP Distribution}
\description{Computes the density, cumulative probability, quantiles, and
random draws for the Zero-Inflated COM-Poisson distribution.}

\usage{
dzicmp(x, lambda, nu, p, z = NULL, max = 100, log = FALSE) 
pzicmp(x, lambda, nu, p, max = 100) 
qzicmp(q, lambda, nu, p, max = 100, log.p = FALSE) 
rzicmp(n, lambda, nu, p, max = 100)
}

\arguments{
  \item{x}{vector of quantiles.}
  \item{q}{vector of probabilities.}
  \item{n}{number of observations.}
  \item{z}{normalizing constant. Can be passed in to save computation;
           otherwise computed internally.}
  \item{lambda}{rate parameter.}
  \item{nu}{dispersion parameter.}
  \item{p}{zero-inflation probability parameter.}
  \item{max}{maximum number to use for truncating infinite sums.}
  \item{log, log.p}{logical; if TRUE, probabilities p are given as log(p).}
}

\value{
\code{dzicmp} gives the density,
\code{pzicmp} gives the cumulative probability,
\code{qzicmp} gives the quantile value, and
\code{rzicmp} generates random numbers.
}

\references{
Kimberly F. Sellers and Andrew~M. Raim (2016). A Flexible Zero-Inflated Model
to Address Data Dispersion. Computational Statistics and Data Analysis, 99,
68-80.
}
\author{Kimberly Sellers, Andrew Raim}

\keyword{ZICMP distribution}
