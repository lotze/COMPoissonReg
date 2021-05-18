#' ZICMP Distribution
#' 
#' Computes the density, cumulative probability, quantiles, and
#' random draws for the zero-inflated COM-Poisson distribution.
#' 
#' @param x vector of quantiles.
#' @param q vector of probabilities.
#' @param n number of observations.
#' @param lambda rate parameter.
#' @param nu dispersion parameter.
#' @param p zero-inflation probability parameter.
#' @param log logical; if TRUE, probabilities are returned on log-scale.
#' @param log.p logical; if TRUE, probabilities p are given as \eqn{\log(p)}.
#' @param method a string: \code{hybrid}, \code{approx}, or \code{trunc}.
#' 
#' @return
#' \describe{
#' \item{dzicmp}{gives the density,}
#' \item{pzicmp}{gives the cumulative probability,}
#' \item{qzicmp}{gives the quantile value, and}
#' \item{rzicmp}{generates random numbers,}
#' \item{ezicmp}{gives the expected value.}
#' \item{vzicmp}{gives the variance.}
#' }
#' 
#' @references
#' Kimberly F. Sellers and Andrew M. Raim (2016). A Flexible Zero-Inflated Model
#' to Address Data Dispersion. Computational Statistics and Data Analysis, 99,
#' 68-80.
#' 
#' @author Kimberly Sellers, Andrew Raim
#' @name ZICMP Distribution
NULL

#' @name ZICMP Distribution
#' @export
dzicmp = function(x, lambda, nu, p, log = FALSE)
{
	n = length(x)
	prep = prep.zicmp(n, lambda, nu, p)
	fx = prep$p*(x==0) + (1-prep$p)*dcmp(x, prep$lambda, prep$nu)
	if (log) { return(log(fx)) } else { return(fx) }
}

#' @name ZICMP Distribution
#' @export
rzicmp = function(n, lambda, nu, p)
{
	prep = prep.zicmp(n, lambda, nu, p)
	s = rbinom(n, size = 1, prob = prep$p)
	(1-s) * rcmp(n, prep$lambda, prep$nu)
}

#' @name ZICMP Distribution
#' @export
pzicmp = function(x, lambda, nu, p)
{
	n = length(x)
	prep = prep.zicmp(n, lambda, nu, p)
	prep$p*(x >= 0) + (1-prep$p)*pcmp(x, prep$lambda, prep$nu)
}

#' @name ZICMP Distribution
#' @export
qzicmp = function(q, lambda, nu, p, log.p = FALSE)
{
	n = length(q)
	prep = prep.zicmp(length(q), lambda, nu, p)

	ymax = getOption("COMPoissonReg.ymax")
	hybrid.tol = getOption("COMPoissonReg.hybrid.tol")
	truncate.tol = getOption("COMPoissonReg.truncate.tol")

	if (log.p) { lq = q } else { lq = log(q) }

	# As in qcmp, check if any requested quantiles are too close to 1. Note
	# that the criteria here depends on p, so we do not include it in the
	# message.
	idx_warn = which(lq > log((1 - truncate.tol)*(1 - prep$p) + prep$p))
	if (length(idx_warn) > 0) {
		msg = sprintf(paste(
			"At least one requested quantile was very close to 1 In",
			"particular, %d of the given probabilities were greater than",
			"(1 - truncate.tol) * (1-p) + p, where truncate_tol = %g.",
			"Associated results may be a consequence of truncation and not",
			"actual quantiles. Consider changing the options",
			"COMPoissonReg.ymax and COMPoissonReg.truncate.tol or reducing",
			"logq"),
			length(idx_warn), truncate.tol)
		warning(msg)
	}

	q_zicmp(lq, prep$lambda, prep$nu, prep$p, hybrid_tol = hybrid.tol,
		truncate_tol = truncate.tol, ymax = ymax)
}

#' @name ZICMP Distribution
#' @export
ezicmp = function(lambda, nu, p, method = "hybrid")
{
	(1-p) * ecmp(lambda, nu, method)
}

#' @name ZICMP Distribution
#' @export
vzicmp = function(lambda, nu, p, method = "hybrid")
{
	ee  = ecmp(lambda, nu, method)
	vv  = vcmp(lambda, nu, method)
	(1-p) * (p*ee^2 + vv)
}

# Extend lambda, nu, and p vectors to be compatible lengths.
# If all are length 1, do not extend them - this is a special case which
# is handled more efficiently.
#
# This also seems to be a good place to make sure parameters are in the right
# space. TBD: do we need to check in the underlying C++ functions?
prep.zicmp = function(n, lambda, nu, p = 0)
{
	L = max(length(lambda), length(nu), length(p))

	stopifnot(all(lambda > 0))
	stopifnot(all(nu > 0))
	stopifnot(all(p >= 0 & p <= 1))

	if (n > 1 && L > 1) { stopifnot(n == L) }
	if (length(lambda) == 1 && L > 1) { lambda = rep(lambda, L) }
	if (length(nu) == 1 && L > 1) { nu = rep(nu, L) }
	if (length(p) == 1 && L > 1) { p = rep(p, L) }
	if (L > 1) { type = "indep" } else { type = "iid" }

	list(lambda = lambda, nu = nu, p = p, type = type)
}
