\name{glm.cmp}
\alias{AIC.cmp}
\alias{BIC.cmp}
\alias{coef.cmp}
\alias{deviance.cmp}
\alias{equitest.cmp}
\alias{glm.cmp}
\alias{leverage.cmp}
\alias{logLik.cmp}
\alias{nu.cmp}
\alias{parametric_bootstrap.cmp}
\alias{predict.cmp}
\alias{print.cmp}
\alias{residuals.cmp}
\alias{sdev.cmp}
\alias{summary.cmp}

\title{Functions to estimate COM-Poisson model parameters (betas and nu)}
\description{ Estimates maximum likelihood estimates for betas and nu using nlminb. }

\usage{
	glm.cmp(formula, initial.est = NULL, nuinit = 1, max = 100, ...)

	\method{AIC}{cmp}(object, ..., k = 2)
	\method{BIC}{cmp}(object, ...)
	\method{coef}{cmp}(object, ...)
	\method{deviance}{cmp}(object, ...)
	\method{equitest}{cmp}(object, ...)
	\method{leverage}{cmp}(object, ...)
	\method{logLik}{cmp}(object, ...)
	\method{nu}{cmp}(object, ...)
	\method{parametric_bootstrap}{cmp}(object, reps = 1000, report.period = reps + 1, ...)
	\method{predict}{cmp}(object, newdata = NULL, ...)
	\method{print}{cmp}(x, ...)
	\method{residuals}{cmp}(object, type = c("raw", "quantile"), ...)
	\method{sdev}{cmp}(object, ...)
	\method{summary}{cmp}(object, ...)
}

\arguments{
  \item{formula}{formula for the COM-Poisson model}
  \item{initial.est}{initial vector of betas, b0_1, ..., b0_p; if NULL, estimated using Poisson GLM}
  \item{nuinit}{initial value for dispersion parameter}
  \item{max}{maximum number to use for truncating infinite sums}
  \item{...}{other model parameters, such as data}
  \item{object}{object of type 'cmp'}
  \item{x}{object of type 'cmp'}
  \item{k}{Penalty per parameter to be used in AIC calculation.}
  \item{newdata}{New covariates to be used for prediction.}
  \item{type}{Type of residual to be computed.}
  \item{reps}{Number of bootstrap repetitions.}
  \item{report.period}{Report progress every \code{report.period} iterations.}
}
\details{
  \code{glm.cmp} finds the maximum likelihood estimates for COM-Poisson.
}

\value{
 An object of class "cmp", from which the coefficients and other information can be computed.
}

\references{
Kimberly F. Sellers & Galit Shmueli (2010). A Flexible Regression Model for
Count Data. Annals of Applied Statistics, 4(2), 943-961.
}
\author{Kimberly Sellers}

\keyword{maximum likelihood estimation}
\keyword{regression}
\keyword{estimate}
