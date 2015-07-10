\name{CMP.MSE}
\alias{CMP.MSE}
\title{COM-Poisson Mean-squared Error}
\description{Computes the Mean-Squared Error associated with the estimated residuals from the COM-Poisson regression.}
\usage{CMP.MSE(CMPbetas, CMPnu, x, y)}

\arguments{
\item{CMPbetas}{Vector estimate of beta-vector, as determined via the COM-Poisson distribution}
\item{CMPnu}{Dispersion estimate determined via the COM-Poisson distribution}
\item{x}{matrix of size nxp}
\item{y}{column vector of length n}
}

\details{ CMP.MSE computes the mean-squared error associated with
the estimated residuals determined from a COM-Poisson regression relating
the matrix x and vector y.}

\value{numeric, mean-squared error associated with
the estimated residuals determined from a COM-Poisson regression relating
the matrix x and vector y.}

\references{A Flexible Regression Model for Count Data, by Sellers & Shmueli, http://ssrn.com/abstract=1127359}
\author{Kimberly Sellers}

\seealso{CMPStdErrors}

\keyword{regression}
\keyword{standard error}
\keyword{mean-squared error}
