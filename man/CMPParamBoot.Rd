\name{CMPParamBoot}
\alias{CMPParamBoot}

\title{Parametric COM-Poisson Bootstrap}
\description{
  Generates n sample estimates for beta and nu.
}
\usage{ 
CMPParamBoot(x, poissonest, betahat, nuhat, n=1000, report.period = n+1,
	max = 100) }

\arguments{
  \item{x}{explanatory variable data}
  \item{poissonest}{initial estimates from Poisson regression}
  \item{betahat}{COM-Poisson estimates from full data set}
  \item{nuhat}{COM-Poisson dispersion estimate from full data set}
  \item{n}{number of sample estimates to create}
  \item{report.period}{Report progress every \code{report.period} iterations}
  \item{max}{maximum number to use for truncating infinite sums}
}

\details{ CMPParamBoot generates n columns of simulated output
(y) data, and uses those simulations along with the explanatory (x)
data to produce n estimates for beta and nu. }

\value{
  \item{Ystar}{a matrix of values with n columns}
  \item{CMPresult}{a matrix of beta and nu estimates containing n rows}
}

\references{
Kimberly F. Sellers & Galit Shmueli (2010). A Flexible Regression Model for
Count Data. Annals of Applied Statistics, 4(2), 943-961.
}
\author{Kimberly Sellers}


\seealso{makeCMPdata, ComputeBetasAndNuHat}

\keyword{COM-Poisson Parametric Bootstrap}
