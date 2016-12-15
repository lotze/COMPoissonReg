\name{glm.zicmp}
\alias{AIC.zicmp}
\alias{BIC.zicmp}
\alias{coef.zicmp}
\alias{deviance.zicmp}
\alias{equitest.zicmp}
\alias{glm.zicmp}
\alias{leverage.zicmp}
\alias{logLik.zicmp}
\alias{nu.zicmp}
\alias{parametric_bootstrap.zicmp}
\alias{predict.zicmp}
\alias{print.zicmp}
\alias{residuals.zicmp}
\alias{sdev.zicmp}
\alias{summary.zicmp}


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
	\method{parametric_bootstrap}{zicmp}(object, reps = 1000, report.period = reps + 1, ...)
	\method{predict}{zicmp}(object, newdata = NULL, ...)
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
  \item{object}{object of type 'zicmp'.}
  \item{x}{object of type 'zicmp'.}
  \item{k}{Penalty per parameter to be used in AIC calculation.}
  \item{newdata}{New covariates to be used for prediction.}
  \item{type}{Type of residual to be computed.}
  \item{reps}{Number of bootstrap reps.}
  \item{report.period}{Report progress every \code{report.period} reps.}
}

\value{
\code{glm.zicmp} produces an object of class "zicmp", from which
coefficients and other information can be computed.
}

\references{
Kimberly F. Sellers & Galit Shmueli (2010). A Flexible Regression Model for
Count Data. Annals of Applied Statistics, 4(2), 943-961.

Kimberly F. Sellers and Andrew~M. Raim (2016). A Flexible Zero-Inflated Model
to Address Data Dispersion. Computational Statistics and Data Analysis, 99,
68-80.
}
\author{Kimberly Sellers}

 \seealso{
       \code{\link{COMPoissonReg}}
     }

\keyword{maximum likelihood estimation}
\keyword{regression}
\keyword{estimate}
