\name{constantCMPfitsandresids}
\alias{constantCMPfitsandresids}

\title{Determining COM-Poisson Fitted Values and Residuals}

\description{
  Computes estimated lambda values, fitted values, and residuals.
} 

\usage{ constantCMPfitsandresids(betahat, nuhat, x, y) }

\arguments{
  \item{betahat}{Maximum likelihood estimate for beta}
  \item{nuhat}{Maximum likelihood estimate for nu}
  \item{x}{matrix of size nxp}
  \item{y}{column vector of length n}
} 

\details{
Computes estimated fitted values, and associated residuals from a COM-Poisson regression with constant dispersion parameter, nu. n=number of records and p=number of predictors.
} 

\value{
   A numeric array of suitable size.
} 

\references{
Kimberly F. Sellers & Galit Shmueli (2010). A Flexible Regression Model for
Count Data. Annals of Applied Statistics, 4(2), 943-961.
}
\author{Kimberly Sellers}

\examples{ # Compute Poisson estimates
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
\keyword{residuals}
\keyword{fitted values}
