\name{CMPStdErrors} \alias{CMPStdErrors}

\title{Standard Errors Associated With COM-Poisson Estimates}

\description{
  This program computes the standard errors associated with the COM-Poisson estimates for the regression and dispersion parameters.
}

\usage{ CMPStdErrors(x, beta, nu, max = 100) }

\arguments{
  \item{x}{matrix of size nxp}
  \item{beta}{Vector estimate determined via the COM-Poisson regression}
  \item{nu}{Dispersion estimate determined via the COM-Poisson regression}
  \item{max}{maximum number to use for truncating infinite sums}
}

\details{
   CMPStdErrors computes the standard errors associated with the COM-Poisson estimates.
}

\value{
  Standard errors.
}

\references{A Flexible Regression Model for Count Data, by Sellers & Shmueli, http://ssrn.com/abstract=1127359}
\author{Kimberly Sellers}

\seealso{CMP.MSE}

\keyword{regression}
\keyword{standard error}
