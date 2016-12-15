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

\title{Functions to fit Zero-Inflated COM-Poisson model}
\description{ Compute maximum likelihood estimates for parameters using optim. }

\usage{
	glm.zicmp(formula.lambda, formula.nu = NULL, formula.p = NULL,
        beta.init = NULL, gamma.init = NULL, zeta.init = NULL, max = 100, ...)

	\method{AIC}{zicmp}(object, ..., k = 2)
	\method{BIC}{zicmp}(object, ...)
	\method{coef}{zicmp}(object, ...)
	\method{deviance}{zicmp}(object, ...)
	\method{equitest}{zicmp}(object, ...)
	\method{leverage}{zicmp}(object, ...)
	\method{logLik}{zicmp}(object, ...)
	\method{nu}{zicmp}(object, ...)
	\method{parametric_bootstrap}{zicmp}(object, ...)
	\method{predict}{zicmp}(object, ...)
	\method{print}{zicmp}(x, ...)
	\method{residuals}{zicmp}(object, type = c("raw", "quantile"), ...)
	\method{sdev}{zicmp}(object, ...)
	\method{summary}{zicmp}(object, ...)
}

\arguments{
  \item{formula.lambda}{regression formula linked to \code{log(lambda)}}
  \item{formula.nu}{regression formula linked to \code{log(nu)}. If NULL, is taken to be intercept only.}
  \item{formula.p}{regression formula linked to \code{logit(p)}. If NULL, is taken to be intercept only.}
  \item{beta.init}{initial values for regression coefficients of \code{lambda}.}
  \item{gamma.init}{initial values for regression coefficients of \code{nu}.}
  \item{zeta.init}{initial values for regression coefficients of \code{p}.}
  \item{max}{maximum number to use for truncating infinite sums.}
  \item{...}{other model parameters, such as data}
  \item{object}{object of type 'zicmp'}
  \item{x}{object of type 'zicmp'}
}

\value{
\code{glm.zicmp} produces an object of class "zicmp", from which
coefficients and other information can be computed.
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
