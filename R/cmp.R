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
#' @param method a string: \code{hybrid}, \code{approx}, or \code{trunc}.
#' 
#' @return 
#' \describe{
#' \item{dcmp}{gives the density,}
#' \item{pcmp}{gives the cumulative probability,}
#' \item{qcmp}{gives the quantile function,}
#' \item{rcmp}{generates random values,}
#' \item{ecmp}{gives the expected value,}
#' \item{vcmp}{gives the variance,}
#' \item{ncmp}{gives the value of the normalizing constant, and}
#' \item{tcmp}{gives the upper value to which we would truncate (see details).}
#' }
#' 
#' @details 
#' TBD: tcmp
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
	prep = prep.zicmp(n, lambda, nu)

	ymax = getOption("COMPoissonReg.ymax")
	hybrid.tol = getOption("COMPoissonReg.hybrid.tol")
	truncate.tol = getOption("COMPoissonReg.truncate.tol")

	if (prep$type == "iid") {
		# Independent and identically distributed case
		out = d_cmp(x, prep$lambda, prep$nu, take_log = log, normalize = TRUE,
			hybrid_tol = hybrid.tol, truncate_tol = truncate.tol, ymax = ymax)
	} else {
		# Independent but not identically distributed case
		out = numeric(n)
		for (i in 1:n) {
			out[i] = d_cmp(x[i], prep$lambda[i], prep$nu[i], take_log = log,
				normalize = TRUE, hybrid_tol = hybrid.tol,
				truncate_tol = truncate.tol, ymax = ymax)
		}
	}

	return(out)
}

#' @name CMP Distribution
#' @export
rcmp = function(n, lambda, nu)
{
	prep = prep.zicmp(n, lambda, nu)

	ymax = getOption("COMPoissonReg.ymax")
	hybrid.tol = getOption("COMPoissonReg.hybrid.tol")
	truncate.tol = getOption("COMPoissonReg.truncate.tol")

	if (prep$type == "iid") {
		# Independent and identically distributed case
		out = r_cmp(n, prep$lambda, prep$nu, hybrid.tol, truncate.tol, ymax)
	} else {
		# Independent but not identically distributed case
		out = numeric(n)
		for (i in 1:n) {
			out[i] = r_cmp(1, prep$lambda[i], prep$nu[i], hybrid.tol, truncate.tol, ymax)
		}
	}

	return(out)
}

#' @name CMP Distribution
#' @export
pcmp = function(x, lambda, nu)
{
	n = length(x)
	prep = prep.zicmp(n, lambda, nu)

	ymax = getOption("COMPoissonReg.ymax")
	hybrid.tol = getOption("COMPoissonReg.hybrid.tol")
	truncate.tol = getOption("COMPoissonReg.truncate.tol")

	if (prep$type == "iid") {
		# Independent and identically distributed case
		out = p_cmp(x, prep$lambda, prep$nu, hybrid_tol = hybrid.tol,
			truncate_tol = truncate.tol, ymax = ymax)
	} else {
		# Independent but not identically distributed case
		out = numeric(n)
		for (i in 1:n) {
			out[i] = p_cmp(x[i], prep$lambda[i], prep$nu[i],
				hybrid_tol = hybrid.tol, truncate_tol = truncate.tol,
				ymax = ymax)
		}
	}

	return(out)
}

#' @name CMP Distribution
#' @export
qcmp = function(q, lambda, nu, log.p = FALSE)
{
	n = length(q)
	prep = prep.zicmp(n, lambda, nu)

	ymax = getOption("COMPoissonReg.ymax")
	hybrid.tol = getOption("COMPoissonReg.hybrid.tol")
	truncate.tol = getOption("COMPoissonReg.truncate.tol")

	if (log.p) { lq = q } else { lq = log(q) }

	if (prep$type == "iid") {
		# Independent and identically distributed case
		out = q_cmp(lq, prep$lambda, prep$nu, hybrid_tol = hybrid.tol,
			truncate_tol = truncate.tol, ymax = ymax)
	} else {
		# Independent but not identically distributed case
		out = numeric(n)
		for (i in 1:n) {
			out[i] = q_cmp(lq[i], prep$lambda[i], prep$nu[i],
				hybrid_tol = hybrid.tol, truncate_tol = truncate.tol, ymax = ymax)
		}
	}

	return(out)
}

#' @name CMP Distribution
#' @export
ecmp = function(lambda, nu, method = "hybrid")
{
	# If lambda and nu are vectors, assume they are not repeats of the same
	# value and do not attempt to save time if this is the case.

	n = max(length(lambda), length(nu))
	prep = prep.zicmp(n, lambda, nu)

	ymax = getOption("COMPoissonReg.ymax")
	hybrid.tol = getOption("COMPoissonReg.hybrid.tol")
	truncate.tol = getOption("COMPoissonReg.truncate.tol")

	# For hybrid method, use asymptotic approximation if this is true.
	# Otherwise use truncated sum.
	is.approx.valid = -1/prep$nu * log(prep$lambda) < log(hybrid.tol)

	out = numeric(n)
	for (i in 1:n) {
		if (method == "approx" || (method == "hybrid" && is.approx.valid[i])) {
			# Use 1st derivative of the log-normalizing constant to compute the
			# expected value.
			out[i] = prep$lambda[i] * grad(ncmp, prep$lambda[i],
				nu = prep$nu[i], log = TRUE)
		} else if (method == "trunc" || (method == "hybrid" && !is.approx.valid[i])) {
			# Compute the expected value by a simple truncated sum
			x.seq = seq(0, tcmp(prep$lambda[i], prep$nu[i]))
			out[i] = sum(x.seq * dcmp(x.seq, prep$lambda[i], prep$nu[i]))
		} else {
			msg = sprintf("Method not recognized: %s", method)
			stop(msg)
		}
	}


	return(out)
}

#' @name CMP Distribution
#' @export
vcmp = function(lambda, nu, method = "hybrid")
{
	# If lambda and nu are vectors, assume they are not repeats of the same
	# value and do not attempt to save time if this is the case.
	n = max(length(lambda), length(nu))
	prep = prep.zicmp(n, lambda, nu)

	ymax = getOption("COMPoissonReg.ymax")
	hybrid.tol = getOption("COMPoissonReg.hybrid.tol")
	truncate.tol = getOption("COMPoissonReg.truncate.tol")

	# For hybrid method, use asymptotic approximation if this is true.
	# Otherwise use truncated sum.
	is.approx.valid = -1/prep$nu * log(prep$lambda) < log(hybrid.tol)

	# Compute expected value
	ev = ecmp(lambda, nu, method)

	out = numeric(n)
	for (i in 1:n) {
		if (method == "approx" || (method == "hybrid" && is.approx.valid[i])) {
			# Use 2nd derivative of the log-normalizing constant to compute the
			# variance. The expression for dd is equal to Var(X) - E(X). Note
			# that we use `method = hybrid` for ncmp itself. This can be changed
			# via the global options, if desired.
			dd = prep$lambda[i]^2 * hessian(ncmp, prep$lambda[i],
			 	nu = prep$nu[i], log = TRUE)
			out[i] = dd + ev[i]
		} else if (method == "trunc" || (method == "hybrid" && !is.approx.valid[i])) {
			# Compute the expected value by a simple truncated sum
			x.seq = seq(0, tcmp(prep$lambda[i], prep$nu[i]))
			out[i] = sum(x.seq^2 * dcmp(x.seq, prep$lambda[i], prep$nu[i])) - ev[i]^2
		} else {
			msg = sprintf("Method not recognized: %s", method)
			stop(msg)
		}
	}

	return(out)
}

#' @name CMP Distribution
#' @export
ncmp = function(lambda, nu, log = FALSE, method = "hybrid")
{
	# If lambda and nu are vectors, assume they are not repeats of the same
	# value and do not attempt to save time if this is the case.
	n = max(length(lambda), length(nu))
	prep = prep.zicmp(n, lambda, nu)

	ymax = getOption("COMPoissonReg.ymax")
	hybrid.tol = getOption("COMPoissonReg.hybrid.tol")
	truncate.tol = getOption("COMPoissonReg.truncate.tol")

	if (method == "trunc") {
		out = z_trunc(prep$lambda, prep$nu, tol = truncate.tol, take_log = log,
			ymax = ymax)
	} else if (method == "approx") {
		out = z_approx(prep$lambda, prep$nu, take_log = log)
	} else if (method == "hybrid") {
		out = z_hybrid(prep$lambda, prep$nu, take_log = log,
			hybrid_tol = hybrid.tol, truncate_tol = truncate.tol, ymax = ymax)
	} else {
		msg = sprintf("Method not recognized: %s", method)
		stop(msg)
	}

	return(out)
}

#' @name CMP Distribution
#' @export
tcmp = function(lambda, nu)
{
	ymax = getOption("COMPoissonReg.ymax")
	truncate.tol = getOption("COMPoissonReg.truncate.tol")

	n = max(length(lambda), length(nu))
	prep = prep.zicmp(n, lambda, nu)
	y_trunc(prep$lambda, prep$nu, tol = truncate.tol, ymax = ymax)

}
