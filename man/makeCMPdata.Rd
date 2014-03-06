\name{makeCMPdata}
\alias{makeCMPdata}

\title{Create Simulated COM-Poisson Data}
\description{ Creates simulated data from a COM-Poisson distribution. } 

\usage{ makeCMPdata(x, beta, nu) }

\arguments{
  \item{x}{matrix of explanatory variables}
  \item{beta}{parameter estimates associated with x}
  \item{nu}{estimated dispersion parameter}
} 

\details{
 makeCMPdata produces simulated output data (y) from a COM-Poisson distribution with estimated parameters beta and nu.
} 

\value{
 Vector of simulated y values.
} 

\references{A Flexible Regression Model for Count Data, by Sellers & Shmueli, http://ssrn.com/abstract=1127359}
\author{Kimberly Sellers}

\seealso{CMPParamBoot, ComputeBetasAndNuHat} 

\examples{
}

\keyword{COM-Poisson simulation}
