\name{CMP Distribution}
\alias{rcmp}
\alias{dcmp}
\alias{pcmp}
\alias{qcmp}

\title{COM-Poisson Distribution}
\description{Functions for the COM-Poisson distribution.}

\usage{
dcmp(x, lambda, nu, z = NULL, log = FALSE, max = 100)
pcmp(x, lambda, nu, max = 100)
qcmp(q, lambda, nu, max = 100, log.p = FALSE)
rcmp(n, lambda, nu, max = 100)
}

\arguments{
  \item{x}{vector of quantiles.}
  \item{q}{vector of probabilities.}
  \item{n}{number of observations.}
  \item{z}{normalizing constant. Can be passed in to save computation;
           otherwise computed internally.}
  \item{lambda}{rate parameter.}
  \item{nu}{dispersion parameter.}
  \item{max}{maximum number to use for truncating infinite sums.}
  \item{log, log.p}{logical; if TRUE, probabilities p are given as log(p).}
}

\value{
\code{dcmp} gives the density,
\code{pcmp} gives the cumulative distribution function,
\code{qcmp} gives the quantile function, and
\code{rcmp} generates random deviates.
}

\references{
Kimberly F. Sellers & Galit Shmueli (2010). A Flexible Regression Model for
Count Data. Annals of Applied Statistics, 4(2), 943-961.
}
\author{Kimberly Sellers}

\keyword{COM-Poisson distribution}
