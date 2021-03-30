#' COM-Poisson Distribution
#' 
#' Functions for the COM-Poisson distribution.
#' 
#' @param x vector of quantiles.
#' @param q vector of probabilities.
#' @param n number of observations.
#' @param lambda rate parameter.
#' @param nu dispersion parameter.
#' @param log logical; if TRUE, probabilities are returned on log-scale.
#' @param log.p logical; if TRUE, probabilities \code{p} are given as \eqn{\log(p)}.
#' 
#' @return 
#' \describe{
#' \item{dcmp}{gives the density,}
#' \item{pcmp}{gives the cumulative probability,}
#' \item{qcmp}{gives the quantile function,}
#' \item{rcmp}{generates random values, and}
#' \item{ecmp}{gives the expected value,}
#' \item{vcmp}{gives the variance.}
#' }
#' 
#' @references
#' Kimberly F. Sellers & Galit Shmueli (2010). A Flexible Regression Model for
#' Count Data. Annals of Applied Statistics, 4(2), 943-961.
#' 
#' @author Kimberly Sellers
#' @keywords COM-Poisson distribution
#' @name CMP Distribution
NULL

#' @name CMP Distribution
#' @export
dcmp = function(x, lambda, nu, log = FALSE)
{
	prep = prep.zicmp(length(x), lambda, nu)
	dcmp_cpp(x, prep$lambda, prep$nu, take_log = log)
}

#' @name CMP Distribution
#' @export
rcmp = function(n, lambda, nu)
{
	prep = prep.zicmp(n, lambda, nu)
	ymax = getOption("COMPoissonReg.ymax", default = 1e6)
	rcmp_cpp(n, prep$lambda, prep$nu, ymax = ymax)
}

#' @name CMP Distribution
#' @export
pcmp = function(x, lambda, nu)
{
	prep = prep.zicmp(length(x), lambda, nu)
	pcmp_cpp(x, prep$lambda, prep$nu)
}

#' @name CMP Distribution
#' @export
qcmp = function(q, lambda, nu, log.p = FALSE)
{
	prep = prep.zicmp(length(q), lambda, nu)
	if (log.p) {
		log.q = q
	} else {
		log.q = log(q)
	}
	ymax = getOption("COMPoissonReg.ymax", default = 1e6)
	qcmp_cpp(log.q, prep$lambda, prep$nu, ymax = ymax)
}

#' @name CMP Distribution
#' @export
ecmp = function(lambda, nu)
{
	n = max(length(lambda), length(nu))
	prep = prep.zicmp(n, lambda, nu)

	res = numeric(n)
	for (i in 1:n) {
		res[i] = prep$lambda[i] * grad.fwd(z_hybrid, prep$lambda[i], nu = prep$nu[i], take_log = TRUE)
	}
	return(res)
}

#' @name CMP Distribution
#' @export
vcmp = function(lambda, nu)
{
	n = max(length(lambda), length(nu))
	ymax = getOption("COMPoissonReg.ymax", default = 1e6)
	prep = prep.zicmp(n, lambda, nu)

	# These tolerances will be made into options in a later version of the
	# code. For now, too many places need to be changed from using default
	# values.
	hybrid_tol = 1e-2
	truncate_tol = 1e-6

	ev = ecmp(lambda, nu)
	if (length(ev) == 1) { ev = rep(ev, n) }

	out = numeric(n)
	for (i in 1:n) {
		# This expression for dd is equal to Var(X) - E(X)
		dd = prep$lambda[i]^2 * hess.fwd(z_hybrid, prep$lambda[i],
		 	nu = prep$nu[i], take_log = TRUE, tol1 = hybrid_tol,
		 	tol2 = truncate_tol)
		out[i] = dd + ev[i]
	}

	return(out)
}
