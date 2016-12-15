\name{glm.cmp}
\alias{glm.cmp}
\alias{deviance.cmp}
\alias{equitest.cmp}
\alias{nu.cmp}
\alias{sdev.cmp}
\alias{leverage.cmp}
\alias{parametric_bootstrap.cmp}
\alias{predict.cmp}
\alias{residuals.cmp}

\title{Functions to estimate COM-Poisson model parameters (betas and nu)}
\description{ Estimates maximum likelihood estimates for betas and nu using nlminb. }

\usage{
	glm.cmp(formula, initial.est=NULL, nuinit=1, max=100, ...)

	## S3 function for class 'cmp'
	\method{nu}{cmp}(object, ...)

	## S3 function for class 'cmp'
	\method{equitest}{cmp}(object, ...)

	## S3 function for class 'cmp'
	\method{sdev}{cmp}(object, ...)

	## S3 function for class 'cmp'
	\method{deviance}{cmp}(object, ...)

	## S3 function for class 'cmp'
	\method{leverage}{cmp}(object, ...)

	## S3 function for class 'cmp'
	\method{predict}{cmp}(object, ...)

	## S3 function for class 'cmp'
	\method{residuals}{cmp}(object, ...)

	## S3 function for class 'cmp'
	\method{parametric_bootstrap}{cmp}(object, ...)
}

\arguments{
  \item{formula}{formula for the COM-Poisson model}
  \item{initial.est}{initial vector of betas, b0_1, ..., b0_p; if NULL, estimated using Poisson GLM}
  \item{nuinit}{initial value for nu}
  \item{max}{maximum number to use for truncating infinite sums}
  \item{...}{other model parameters, such as data}
  \item{object}{object of type 'cmp'}
}
\details{
  \code{glm.cmp} finds the maximum likelihood estimates for the COM-Poisson estimators (betas and nu) using the function ComputeBetasAndNuHat.
}

\value{
 An object of class "cmp", from which the coefficients and other information can be computed.
}

\references{
Kimberly F. Sellers & Galit Shmueli (2010). A Flexible Regression Model for
Count Data. Annals of Applied Statistics, 4(2), 943-961.
}
\author{Kimberly Sellers}

 \seealso{
       \code{\link{COMPoissonReg}}
       \code{\link{ComputeBetasAndNuHat}}
     }

\keyword{maximum likelihood estimation}
\keyword{regression}
\keyword{estimate}
