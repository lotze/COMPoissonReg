\name{CMP Distribution}
\alias{rcmp}
\alias{dcmp}
\alias{pcmp}
\alias{qcmp}

\title{COM-Poisson Distribution}
\description{Functions for the COM-Poisson distribution.}

\usage{
dcmp(x, lambda, nu, log = FALSE)
pcmp(x, lambda, nu)
qcmp(q, lambda, nu, log.p = FALSE)
rcmp(n, lambda, nu)
}

\arguments{
  \item{x}{vector of quantiles.}
  \item{q}{vector of probabilities.}
  \item{n}{number of observations.}
  \item{lambda}{rate parameter.}
  \item{nu}{dispersion parameter.}
  \item{log, log.p}{logical; if TRUE, probabilities p are given as log(p).}
}

\value{
\code{dcmp} gives the density,
\code{pcmp} gives the cumulative probability,
\code{qcmp} gives the quantile function, and
\code{rcmp} generates random values.
}

\references{
Kimberly F. Sellers & Galit Shmueli (2010). A Flexible Regression Model for
Count Data. Annals of Applied Statistics, 4(2), 943-961.
}
\author{Kimberly Sellers}

\keyword{COM-Poisson distribution}
