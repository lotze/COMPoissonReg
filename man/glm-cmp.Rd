\name{glm.cmp}
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
\alias{vcov.zicmp}

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
\alias{vcov.cmp}

\title{COM-Poisson and Zero-Inflated COM-Poisson regression}
\description{Fit COM-Poisson regression using maximum likelihood estimation.
Zero-Inflated COM-Poisson can be fit by specifying a regression for the
overdispersion parameter.

The COM-Poisson regression model is
\deqn{
y_i \sim \rm{CMP}(\lambda_i, \nu_i), \;\;\;
\log \lambda_i = \bm{x}_i^\top \beta, \;\;\;
\log \nu_i = \bm{s}_i^\top \gamma.
}{
y_i ~ CMP(lambda_i, nu_i),
log lambda_i = x_i^T beta,
log nu_i = s_i^T gamma.
}

The Zero-Inflated COM-Poisson regression model assumes that \eqn{y_i} is 0
with probability \eqn{p_i} or \eqn{y_i^*} with probability \eqn{1 - p_i},
where
\deqn{
y_i^* \sim \rm{CMP}(\lambda_i, \nu_i), \;\;\;
\log \lambda_i = \bm{x}_i^\top \beta, \;\;\;
\log \nu_i = \bm{s}_i^\top \gamma, \;\;\;
\log p_i = \bm{w}_i^\top \zeta.
}{
y_i^* ~ CMP(lambda_i, nu_i),
log lambda_i = x_i^T beta,
log nu_i = s_i^T gamma,
log p_i = w_i^T zeta.
}
}

\usage{
	glm.cmp(formula.lambda, formula.nu = ~ 1, formula.p = NULL,
        beta.init = NULL, gamma.init = NULL, zeta.init = NULL, max = 100, ...)

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
	\method{vcov}{cmp}(object, ...)


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
	\method{vcov}{zicmp}(object, ...)
}

\arguments{
  \item{formula.lambda}{regression formula linked to \code{log(lambda)}}
  \item{formula.nu}{regression formula linked to \code{log(nu)}. By default, is taken to be intercept only.}
  \item{formula.p}{regression formula linked to \code{logit(p)}. If NULL (the default), zero-inflation term is excluded from the model.}
  \item{beta.init}{initial values for regression coefficients of \code{lambda}.}
  \item{gamma.init}{initial values for regression coefficients of \code{nu}.}
  \item{zeta.init}{initial values for regression coefficients of \code{p}.}
  \item{max}{maximum number to use for truncating infinite sums.}
  \item{...}{other model parameters, such as data}
  \item{object}{object of type 'cmp' or 'zicmp'.}
  \item{x}{object of type 'cmp' or 'zicmp'.}
  \item{k}{Penalty per parameter to be used in AIC calculation.}
  \item{newdata}{New covariates to be used for prediction.}
  \item{type}{Type of residual to be computed.}
  \item{reps}{Number of bootstrap repetitions.}
  \item{report.period}{Report progress every \code{report.period} iterations.}
}

\value{
\code{glm.cmp} produces an object of either class 'cmp' or 'zicmp', depending
on whether zero-inflation is used in the model. From this object, coefficients
and other information can be extracted.
}

\references{
Kimberly F. Sellers & Galit Shmueli (2010). A Flexible Regression Model for
Count Data. Annals of Applied Statistics, 4(2), 943-961.

Kimberly F. Sellers and Andrew M. Raim (2016). A Flexible Zero-Inflated Model
to Address Data Dispersion. Computational Statistics and Data Analysis, 99,
68-80.
}
\author{Kimberly Sellers, Thomas Lotze, Andrew Raim}

\keyword{maximum likelihood estimation}
\keyword{regression}
\keyword{estimate}
