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
#' \item{ecmp}{gives the expected value.}
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
	n = length(x)
	ymax = getOption("COMPoissonReg.ymax")
	hybrid_tol = getOption("COMPoissonReg.hybrid_tol")
	truncate_tol = getOption("COMPoissonReg.truncate_tol")
	prep = prep.zicmp(n, lambda, nu)

	if (prep$type == "iid") {
		# Independent and identically distributed
		out = d_cmp(x, prep$lambda, prep$nu, take_log = log, normalize = TRUE,
			hybrid_tol = hybrid_tol, truncate_tol = truncate_tol, ymax = ymax)
	} else {
		# Independent but not identically distributed
		x = numeric(n)
		for (i in 1:n) {
			out[i] = d_cmp(x[i], prep$lambda[i], prep$nu[i], take_log = log,
				normalize = TRUE, hybrid_tol = hybrid_tol,
				truncate_tol = truncate_tol, ymax = ymax)
		}
	}

	return(out)
}

#' @name CMP Distribution
#' @export
rcmp = function(n, lambda, nu)
{
	ymax = getOption("COMPoissonReg.ymax")
	hybrid_tol = getOption("COMPoissonReg.hybrid_tol")
	truncate_tol = getOption("COMPoissonReg.truncate_tol")
	prep = prep.zicmp(n, lambda, nu)

	if (prep$type == "iid") {
		# Independent and identically distributed
		out = r_cmp(n, prep$lambda, prep$nu, hybrid_tol, truncate_tol, ymax)
	} else {
		# Independent but not identically distributed
		out = numeric(n)
		for (i in 1:n) {
			out[i] = r_cmp(1, prep$lambda[i], prep$nu[i], hybrid_tol, truncate_tol, ymax)
		}
	}

	return(out)
}

#' @name CMP Distribution
#' @export
pcmp = function(x, lambda, nu)
{
	n = length(x)
	ymax = getOption("COMPoissonReg.ymax")
	hybrid_tol = getOption("COMPoissonReg.hybrid_tol")
	truncate_tol = getOption("COMPoissonReg.truncate_tol")
	prep = prep.zicmp(n, lambda, nu)

	if (prep$type == "iid") {
		# Independent and identically distributed
		out = p_cmp(x, prep$lambda, prep$nu, hybrid_tol, truncate_tol, ymax)
	} else {
		# Independent but not identically distributed
		out = numeric(n)
		for (i in 1:n) {
			out[i] = p_cmp(x[i], prep$lambda[i], prep$nu[i], hybrid_tol, truncate_tol, ymax)
		}
	}

	return(out)
}

#' @name CMP Distribution
#' @export
qcmp = function(q, lambda, nu, log.p = FALSE)
{
	n = length(x)
	ymax = getOption("COMPoissonReg.ymax")
	hybrid_tol = getOption("COMPoissonReg.hybrid_tol")
	truncate_tol = getOption("COMPoissonReg.truncate_tol")
	prep = prep.zicmp(n, lambda, nu)

	if (log.p) { lq = q } else { lq = log(q) }

	if (prep$type == "iid") {
		# Independent and identically distributed
		out = q_cmp(lq, prep$lambda, prep$nu, hybrid_tol, truncate_tol, ymax)
	} else {
		# Independent but not identically distributed
		out = numeric(n)
		for (i in 1:n) {
			out[i] = q_cmp(lq[i], prep$lambda[i], prep$nu[i], hybrid_tol, truncate_tol, ymax)
		}
	}

	return(out)
}

#' @name CMP Distribution
#' @export
ecmp = function(lambda, nu)
{
	n = max(length(lambda), length(nu))
	ymax = getOption("COMPoissonReg.ymax")
	hybrid_tol = getOption("COMPoissonReg.hybrid_tol")
	truncate_tol = getOption("COMPoissonReg.truncate_tol")
	prep = prep.zicmp(n, lambda, nu)

	out = numeric(n)
	for (i in 1:n) {
		out[i] = prep$lambda[i] * grad.fwd(z_hybrid, prep$lambda[i],
			nu = prep$nu[i], take_log = TRUE, hybrid_tol = hybrid_tol,
			truncate_tol = truncate_tol, ymax = ymax)
	}

	return(out)
}

#' @name CMP Distribution
#' @export
vcmp = function(lambda, nu)
{
	n = max(length(lambda), length(nu))
	ymax = getOption("COMPoissonReg.ymax")
	hybrid_tol = getOption("COMPoissonReg.hybrid_tol")
	truncate_tol = getOption("COMPoissonReg.truncate_tol")
	prep = prep.zicmp(n, lambda, nu)

	ev = ecmp(lambda, nu)
	if (length(ev) == 1) { ev = rep(ev, n) }

	out = numeric(n)
	for (i in 1:n) {
		# This expression for dd is equal to Var(X) - E(X)
		dd = prep$lambda[i]^2 * hess.fwd(z_hybrid, prep$lambda[i],
		 	nu = prep$nu[i], take_log = TRUE, hybrid_tol = hybrid_tol,
		 	truncate_tol = truncate_tol, ymax = ymax)
		out[i] = dd + ev[i]
	}

	return(out)
}
