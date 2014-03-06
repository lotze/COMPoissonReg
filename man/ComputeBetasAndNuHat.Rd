\name{ComputeBetasAndNuHat}
\alias{ComputeBetasAndNuHat}
\title{Estimate COM-Poisson parameters (betas and nu)}
\description{ Estimates maximum likelihood estimates for betas and nu using nlminb. } 

\usage{ ComputeBetasAndNuHat(x, y, betainit, nuinit, max) }

\arguments{
  \item{x}{matrix of size nxp}
  \item{y}{column vector of length n}
  \item{betainit}{initial vector of betas, b0_1, ..., b0_p}
  \item{nuinit}{initial value for nu}
  \item{max}{maximum number to use for truncating infinite sums}
} 
\details{
ComputeBetasAndNuHat finds the maximum likelihood estimates for the COM-Poisson estimators (betas and nu) using constrained optimization (nlminb).  The dispersion parameter, nu, is greater than or equal to zero. n=number of records and p=number of predictors.
} 

\value{
 Maximum likelihood estimates for betas and nu using constrained optimization (nlminb).
} 

\references{A Flexible Regression Model for Count Data, by Sellers & Shmueli, http://ssrn.com/abstract=1127359}
\author{Kimberly Sellers}

\examples{
}

\keyword{maximum likelihood estimation}
\keyword{regression}
\keyword{estimate}
