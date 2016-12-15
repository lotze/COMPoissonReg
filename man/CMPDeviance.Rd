\name{CMPDeviance}
\alias{CMPDeviance}

\title{CMP Deviance}
\description{Computes the COM-Poisson deviance.}

\usage{ CMPDeviance(x, y, betahat, nuhat, leverage, max = 150) }

\arguments{
  \item{x}{matrix of size nxp}
  \item{y}{column vector of length n}
  \item{betahat}{COM-Poisson maximum likelihood estimates for beta vector}
  \item{nuhat}{COM-Poisson estimate for dispersion parameter}
  \item{leverage}{COM-Poisson estimate for leverage}
  \item{max}{maximum number used to approximate infinite sum calculation}
}

\details{
  CMPDeviance computes the COM-Poisson deviances exactly, through the use of the associated exact log likelihood computations, and associated leverage computations.
}

\value{
  Exact COM-Poisson deviances.
}

\references{
Kimberly F. Sellers & Galit Shmueli (2010). A Flexible Regression Model for
Count Data. Annals of Applied Statistics, 4(2), 943-961.
}
\author{Kimberly Sellers}

\seealso{leverage}

\keyword{regression}
\keyword{deviance}
