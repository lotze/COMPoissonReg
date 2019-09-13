#' Estimate parameters for COM-Poisson regression
#' 
#' This package offers the ability to compute the parameter estimates
#' for a COM-Poisson or zero-inflated (ZI) COM-Poisson regression and
#' associated standard errors.  This package also provides a hypothesis
#' test for determining statistically significant data dispersion, and
#' other model diagnostics.
#' 
#' @details
#' This package offers the ability to compute COM-Poisson parameter
#' estimates and associated standard errors for a regular regression
#' model or a zero-inflated regression model (via the \code{glm_cmp}
#' function).
#' 
#' Further, the user can perform a hypothesis test to determine the
#' statistically significant need for using COM-Poisson regression
#' to model the data.  The test addresses the matter of statistically
#' significant dispersion.
#' 
#' The main order of functions for COM-Poisson regression is as follows:
#' \enumerate{
#' \item Compute Poisson estimates (using \code{glm} for Poisson regression
#'     or \code{pscl} for ZIP regression).
#' \item Use Poisson estimates as starting values to determine COM-Poisson
#'     estimates (using \code{glm_cmp}).
#' \item Compute associated standard errors (using \code{sdev} function).
#' }
#' 
#' From here, there are many ways to proceed, so order is irrelevant:
#' \itemize{
#' \item Perform a hypothesis test to assess for statistically significant
#'       dispersion (using equitest or parametric_bootstrap).
#' \item Compute leverage (using leverage) and deviance (using deviance).
#' \item Predict the outcome for new examples, using predict.
#' }
#' 
#' The package also supports fitting of the zero-inflated COM-Poisson model
#' (ZICMP). Most of the tools available for COM-Poisson are also available
#' for ZICMP.
#' 
#' As of version 0.5.0 of this package, a hybrid method is used to compute
#' the normalizing constant \eqn{z(\lambda, \nu)} for the COM-Poisson density.
#' A closed-form approximation (Shmueli et al, 2005; Gillispie & Green, 2015)
#' to the exact sum is used if the given \eqn{\lambda} is sufficiently large
#' and \eqn{\nu} is sufficiently small. Otherwise, an exact summation is used,
#' except that the number of terms is truncated to meet a given accuracy.
#' Previous versions of the package used simple truncation (defaulting to 100
#' terms), but this was found to be inaccurate in some settings.
#' 
#' @author Kimberly Sellers, Thomas Lotze, Andrew M. Raim
#' 
#' @references
#' Steven B. Gillispie & Christopher G. Green (2015) Approximating the
#' Conway-Maxwell-Poisson distribution normalization constant, Statistics,
#' 49:5, 1062-1073.
#' 
#' Kimberly F. Sellers & Galit Shmueli (2010). A Flexible Regression Model for
#' Count Data. Annals of Applied Statistics, 4(2), 943-961.
#' 
#' Kimberly F. Sellers and Andrew M. Raim (2016). A Flexible Zero-Inflated
#' Model to Address Data Dispersion. Computational Statistics and Data
#' Analysis, 99, 68-80.
#' 
#' Galit Shmueli, Thomas P. Minka, Joseph B. Kadane, Sharad Borle, and Peter
#' Boatwright (2005). A useful distribution for fitting discrete data: revival
#' of the Conway-Maxwell-Poisson distribution. Journal of Royal Statistical
#' Society C, 54, 127-142.
#' 
#' @examples 
#' ## load freight data
#' data(freight)
#' 
#' # Fit standard Poisson model
#' glm_out = glm(broken ~ transfers, data=freight,
#'   family=poisson, na.action=na.exclude)
#' print(glm_out)
#' 
#' # Fit COM-Poisson model (with intercept-only regression linked to the
#' # dispersion parameter)
#' cmp_out = glm_cmp(broken ~ transfers, data=freight)
#' print(cmp_out)
#' coef(cmp_out)
#' nu(cmp_out)[1]
#' 
#' # Compute associated standard errors
#' sdev(cmp_out)
#' 
#' # Get the full covariance matrix for the estimates
#' vcov(cmp_out)
#' 
#' # Likelihood ratio test for dispersion parameter
#' # Test for H_0: dispersion equal to 1 vs. H_1: not equal to 1
#' # (i.e. Poisson vs. COM-Poisson regression models)
#' lrt = equitest(cmp_out)
#' 
#' # Compute constant COM-Poisson leverage
#' lev = leverage(cmp_out)
#' 
#' \dontrun{
#' # Compute constant COM-Poisson deviances
#' dev = deviance(cmp_out)
#' }
#' 
#' # Compute fitted values
#' y_hat = predict(cmp_out, newdata=freight)
#' 
#' # Compute residual values
#' res = residuals(cmp_out)
#' print(summary(res))
#' 
#' # Compute MSE
#' mean(res^2)
#' 
#' # Compute predictions on new data
#' new_data = data.frame(transfers=(0:10))
#' y_hat = predict(cmp_out, newdata=new_data)
#' plot(0:10, y_hat, type="l",
#'   xlab="number of transfers", ylab="predicted number broken")
#' 
#' \dontrun{
#' # Compute parametric bootstrap results and use them to generate
#' # 0.95 confidence intervals for parameters.
#' cmp_boot = parametric_bootstrap(cmp_out, reps=1000)
#' print(apply(cmp_boot, 2, quantile, c(0.025,0.975)))
#' }
#' 
#' \dontrun{
#' ## load couple data
#' data(couple)
#' 
#' # Fit standard Poisson model
#' glm_out = glm(UPB ~ EDUCATION + ANXIETY, data=couple, family=poisson)
#' print(glm_out)
#' 
#' # Fit ZICMP model
#' zicmp_out = glm_cmp(UPB ~ EDUCATION + ANXIETY,
#'   formula_nu = ~ 1,
#'   formula_p = ~ EDUCATION + ANXIETY,
#'   data=couple)
#' print(zicmp_out)
#' 
#' # Compute standard errors for estimates of coefficients
#' sdev(zicmp_out)
#' 
#' # Get the full covariance matrix for the estimates
#' vcov(zicmp_out)
#' 
#' # Likelihood ratio test for equidispersion (H0: nu = 1 vs H1: not)
#' equitest(zicmp_out)
#' 
#' # Compute fitted values
#' y_hat = predict(zicmp_out)
#' 
#' # Compute residuals
#' res_raw = residuals(zicmp_out, type = "raw")
#' res_quan = residuals(zicmp_out, type = "quantile")
#' print(summary(res_raw))
#' print(summary(res_quan))
#' 
#' # Compute predictions on new data
#' new_data = data.frame(EDUCATION = round(1:20 / 20), ANXIETY = seq(-3,3, length_out = 20))
#' y_hat_new = predict(zicmp_out, newdata=new_data)
#' print(y_hat_new)
#' 
#' # Compute parametric bootstrap results and use them to generate
#' # 0.95 confidence intervals for parameters.
#' zicmp_boot = parametric_bootstrap(zicmp_out, reps=1000)
#' print(apply(zicmp_boot, 2, quantile, c(0.025,0.975)))
#' }
#' 
#' @name COMPoissonReg-package
#' @useDynLib COMPoissonReg, .registration = TRUE
#' @import Rcpp
#' @import stats
#' @docType package
NULL

#' Package options
#' 
#' Global options used by the COMPoissonReg package.
#' 
#' @details
#' \itemize{
#' \item \code{options(COMPoissonReg.optim.method = 'L-BFGS-B')}
#' \item \code{options(COMPoissonReg.optim.control = list(maxit = 150))}
#' \item \code{options(COMPoissonReg.grad.eps = 1e-5)}
#' \item \code{options(COMPoissonReg.hess.eps = 1e-2)}
#' \item \code{options(COMPoissonReg.ymax = 1e6)}
#' }
#' 
#' @param COMPoissonReg.optim.method Optim method to use when computing
#'   maximum likelihood estimates.
#' @param COMPoissonReg.optim.control A list to be passed to \code{control}
#'   when calling \code{optim}. \code{fnscale} will be ignored if specified.
#' @param COMPoissonReg.grad.eps Distance to be used when finite differences
#'   are taken.
#' @param COMPoissonReg.hess.eps Distance to be used when finite second
#'   differences are taken.
#' @param COMPoissonReg.ymax Maximum count value to be considered. Larger
#'   values are truncated.
#' 
#' @name COMPoissonReg-options
NULL
