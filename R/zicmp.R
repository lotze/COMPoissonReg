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
#' @param control a \code{COMPoissonReg.control} object from \code{get.control}
#' or \code{NULL} to use global default.
#' 
#' @return
#' \describe{
#' \item{dzicmp}{density,}
#' \item{pzicmp}{cumulative probability,}
#' \item{qzicmp}{quantiles,}
#' \item{rzicmp}{generate random variates,}
#' \item{ezicmp}{expected value. and}
#' \item{vzicmp}{variance.}
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
dzicmp = function(x, lambda, nu, p, log = FALSE, control = NULL)
{
	n = length(x)
	prep = prep.zicmp(n, lambda, nu, p)
	fx = prep$p*(x==0) + (1-prep$p)*dcmp(x, prep$lambda, prep$nu, control = control)
	if (log) { return(log(fx)) } else { return(fx) }
}

#' @name ZICMP Distribution
#' @export
rzicmp = function(n, lambda, nu, p, control = NULL)
{
	prep = prep.zicmp(n, lambda, nu, p)
	s = rbinom(n, size = 1, prob = prep$p)
	(1-s) * rcmp(n, prep$lambda, prep$nu, control = control)
}

#' @name ZICMP Distribution
#' @export
pzicmp = function(x, lambda, nu, p, control = NULL)
{
	n = length(x)
	prep = prep.zicmp(n, lambda, nu, p)
	prep$p*(x >= 0) + (1-prep$p)*pcmp(x, prep$lambda, prep$nu, control = control)
}

#' @name ZICMP Distribution
#' @export
qzicmp = function(q, lambda, nu, p, log.p = FALSE, control = NULL)
{
	n = length(q)
	prep = prep.zicmp(length(q), lambda, nu, p)

	if (is.null(control)) { control = getOption("COMPoissonReg.control") }
	stopifnot("COMPoissonReg.control" %in% class(control))
	ymax = control$ymax
	hybrid.tol = control$hybrid.tol
	truncate.tol = control$truncate.tol

	if (log.p) { lq = q } else { lq = log(q) }

	# As in qcmp, check if any requested quantiles are too close to 1. Note
	# that the criteria here depends on p, so we do not include it in the
	# message.
	idx_warn = which(lq > log((1 - truncate.tol)*(1 - prep$p) + prep$p))
	if (length(idx_warn) > 0) {
		msg = sprintf(paste(
			"At least one requested quantile was very close to 1. In",
			"particular, %d of the given probabilities were greater than",
			"(1 - truncate.tol) * (1-p) + p, where truncate_tol = %g.",
			"Associated results may be affected by truncation. Consider",
			"adjusting the controls ymax and truncate.tol or reducing logq."),
			length(idx_warn), truncate.tol)
		warning(msg)
	}

	q_zicmp(lq, prep$lambda, prep$nu, prep$p, hybrid_tol = hybrid.tol,
		truncate_tol = truncate.tol, ymax = ymax)
}

#' @name ZICMP Distribution
#' @export
ezicmp = function(lambda, nu, p, control = NULL)
{
	(1-p) * ecmp(lambda, nu, control = control)
}

#' @name ZICMP Distribution
#' @export
vzicmp = function(lambda, nu, p, control = NULL)
{
	ee  = ecmp(lambda, nu, control = control)
	vv  = vcmp(lambda, nu, control = control)
	(1-p) * (p*ee^2 + vv)
}
