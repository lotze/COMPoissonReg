\name{CMPLeverage}
\alias{CMPLeverage}

\title{Computes leverage for each record}
\description{ Computes leverage for each record, for purpose of detecting influential records.}

\usage{ CMPLeverage(x, y, betahat, nuhat, max) }

\arguments{
  \item{x}{matrix of size n by p}
  \item{y}{column vector of length n}
  \item{betahat}{Maximum likelihood estimate for beta vector}
  \item{nuhat}{Maximum likelihood estimate for nu}
  \item{max}{maximum number to use for truncating infinite sums}
}

\details{
  Leverage computes the Hat matrix: H = W^(1/2) X (X'WX)^{-1} X'W^(1/2) , where X is the design matrix and W is the weight matrix defined in
    Shmueli and Sellers (2009). The leverage for record i is then given by the i^(th) diagonal element of H. Note that, like in other non-linear models, H depends not only on the
    covariate values X but also on the estimated parameters (in W).
}

\value{the Hat matrix}


\references{A Flexible Regression Model for Count Data, by Sellers & Shmueli, http://ssrn.com/abstract=1127359}
\author{Kimberly Sellers}

\seealso{ExactDeviance}

\keyword{regression}
