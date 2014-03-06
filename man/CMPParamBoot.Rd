\name{CMPParamBoot}
\alias{CMPParamBoot}

\title{Parametric COM-Poisson Bootstrap}
\description{
  Generates n sample estimates for beta and nu.
} 
\usage{ CMPParamBoot(x, poissonest, betahat, nuhat, n=1000) }

\arguments{
  \item{x}{explanatory variable data}
  \item{poissonest}{initial estimates from Poisson regression}
  \item{betahat}{COM-Poisson estimates from full data set}
  \item{nuhat}{COM-Poisson dispersion estimate from full data set}
  \item{n}{number of sample estimates to create}
} 

\details{ CMPParamBoot generates n columns of simulated output
(y) data, and uses those simulations along with the explanatory (x)
data to produce n estimates for beta and nu. } 

\value{
  \item{Ystar}{a matrix of values with n columns}
  \item{CMPresult}{a matrix of beta and nu estimates containing n rows}
} 

\references{A Flexible Regression Model for Count Data, by Sellers & Shmueli, http://ssrn.com/abstract=1127359}
\author{Kimberly Sellers}


\seealso{makeCMPdata, ComputeBetasAndNuHat} 

\examples{

}

\keyword{COM-Poisson Parametric Bootstrap}
