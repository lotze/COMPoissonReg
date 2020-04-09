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
#' model or a zero-inflated regression model (via the \code{glm.cmp}
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
#'     estimates (using \code{glm.cmp}).
#' \item Compute associated standard errors (using \code{sdev} function).
#' }
#' 
#' From here, there are many ways to proceed, so order is irrelevant:
#' \itemize{
#' \item Perform a hypothesis test to assess for statistically significant
#'       dispersion (using equitest or parametric.bootstrap).
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
#' glm.out = glm(broken ~ transfers, data=freight,
#'   family=poisson, na.action=na.exclude)
#' print(glm.out)
#' 
#' # Fit COM-Poisson model (with intercept-only regression linked to the
#' # dispersion parameter)
#' cmp.out = glm.cmp(broken ~ transfers, data=freight)
#' print(cmp.out)
#' coef(cmp.out)
#' nu(cmp.out)[1]
#' 
#' # Compute associated standard errors
#' sdev(cmp.out)
#' 
#' # Get the full covariance matrix for the estimates
#' vcov(cmp.out)
#' 
#' # Likelihood ratio test for dispersion parameter
#' # Test for H_0: dispersion equal to 1 vs. H_1: not equal to 1
#' # (i.e. Poisson vs. COM-Poisson regression models)
#' lrt = equitest(cmp.out)
#' 
#' # Compute constant COM-Poisson leverage
#' lev = leverage(cmp.out)
#' 
#' \dontrun{
#' # Compute constant COM-Poisson deviances
#' dev = deviance(cmp.out)
#' }
#' 
#' # Compute fitted values
#' y.hat = predict(cmp.out, newdata=freight)
#' 
#' # Compute residual values
#' res = residuals(cmp.out)
#' print(summary(res))
#' 
#' # Compute MSE
#' mean(res^2)
#' 
#' # Compute predictions on new data
#' new.data = data.frame(transfers=(0:10))
#' y.hat = predict(cmp.out, newdata=new.data)
#' plot(0:10, y.hat, type="l",
#'   xlab="number of transfers", ylab="predicted number broken")
#' 
#' \dontrun{
#' # Compute parametric bootstrap results and use them to generate
#' # 0.95 confidence intervals for parameters.
#' cmp.boot = parametric.bootstrap(cmp.out, reps=1000)
#' print(apply(cmp.boot, 2, quantile, c(0.025,0.975)))
#' }
#' 
#' \dontrun{
#' ## load couple data
#' data(couple)
#' 
#' # Fit standard Poisson model
#' glm.out = glm(UPB ~ EDUCATION + ANXIETY, data=couple, family=poisson)
#' print(glm.out)
#' 
#' # Fit ZICMP model
#' zicmp.out = glm.cmp(UPB ~ EDUCATION + ANXIETY,
#'   formula.nu = ~ 1,
#'   formula.p = ~ EDUCATION + ANXIETY,
#'   data=couple)
#' print(zicmp.out)
#' 
#' # Compute standard errors for estimates of coefficients
#' sdev(zicmp.out)
#' 
#' # Get the full covariance matrix for the estimates
#' vcov(zicmp.out)
#' 
#' # Likelihood ratio test for equidispersion (H0: nu = 1 vs H1: not)
#' equitest(zicmp.out)
#' 
#' # Compute fitted values
#' y.hat = predict(zicmp.out)
#' 
#' # Compute residuals
#' res.raw = residuals(zicmp.out, type = "raw")
#' res.quan = residuals(zicmp.out, type = "quantile")
#' print(summary(res.raw))
#' print(summary(res.quan))
#' 
#' # Compute predictions on new data
#' new.data = data.frame(EDUCATION = round(1:20 / 20), ANXIETY = seq(-3,3, length.out = 20))
#' y.hat.new = predict(zicmp.out, newdata=new.data)
#' print(y.hat.new)
#' 
#' # Compute parametric bootstrap results and use them to generate
#' # 0.95 confidence intervals for parameters.
#' zicmp.boot = parametric.bootstrap(zicmp.out, reps=1000)
#' print(apply(zicmp.boot, 2, quantile, c(0.025,0.975)))
#' 
#' # A CMP example with offset terms.
#' cmp.out = glm.cmp(broken ~ transfers + offset(transfers), data=freight)
#' print(cmp.out)
#' coef(cmp.out)
#' nu(cmp.out)[1]
#' 
#' # A ZICMP example with offset terms.
#' zicmp.out = glm.cmp(UPB ~ EDUCATION + ANXIETY + offset(ANXIETY),
#'     formula.nu = ~ offset(ANXIETY),
#'     formula.p = ~ EDUCATION + ANXIETY + offset(ANXIETY),
#'     data=couple)
#' print(zicmp.out)
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