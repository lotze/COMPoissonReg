\name{LRT}
\alias{LRT}

\title{Likelihood ratio test for statistically significant dispersion}
\description{
 Computes the -2loglikelihood value and associated p-value for testing whether Poisson regression is reasonable (H0) versus the COM-Poisson regression (H1).
} 

\usage{ LRT(x, y, betahat0, betahat, nuhat, max) }

\arguments{
  \item{x}{matrix of size nxp}
  \item{y}{column vector of length n}
  \item{betahat0}{Maximum likelihood estimate for beta vector, obtained under the Poisson distribution}
  \item{betahat}{Maximum likelihood estimate for beta vector, obtained under COM-Poisson distribution}
  \item{nuhat}{Maximum likelihood estimate for nu under COM-Poisson distribution}
  \item{max}{maximum number to use for truncating infinite sums}
} 

\details{
  LRT computes the test statistic and p-value for testing H0:nu=1 vs. a two-sided alternative.
    Note that the p-value relies on an asymptotic normal distribution and should therefore be used only
    with sufficiently large n. n=number of records and p=number of predictors.
} 

\value{
 Test statistic and p-value.
} 

\references{A Flexible Regression Model for Count Data, by Sellers & Shmueli, http://ssrn.com/abstract=1127359}
\author{Kimberly Sellers}

\examples{
}

\keyword{regression}
\keyword{likelihood ratio test}
