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

\references{A Flexible Regression Model for Count Data, by Sellers & Shmueli, http://ssrn.com/abstract=1127359}
\author{Kimberly Sellers}

\seealso{leverage}

\examples{	
}

\keyword{regression}
\keyword{deviance}
