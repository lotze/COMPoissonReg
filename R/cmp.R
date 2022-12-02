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
#' @param control a \code{COMPoissonReg.control} object from \code{get.control}
#' or \code{NULL} to use global default.
#' 
#' @return 
#' \describe{
#' \item{dcmp}{density,}
#' \item{pcmp}{cumulative probability,}
#' \item{qcmp}{quantiles,}
#' \item{rcmp}{generate random variates,}
#' \item{ecmp}{expected value,}
#' \item{vcmp}{variance,}
#' \item{ncmp}{value of the normalizing constant, and}
#' \item{tcmp}{upper value used to compute the normalizing constant under
#' truncation method.}
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
dcmp = function(x, lambda, nu, log = FALSE, control = NULL)
{
	n = length(x)
	prep = prep.zicmp(n, lambda, nu)

	if (is.null(control)) {
		control = getOption("COMPoissonReg.control", default = get.control())
	}
	stopifnot("COMPoissonReg.control" %in% class(control))
	ymax = control$ymax
	hybrid.tol = control$hybrid.tol
	truncate.tol = control$truncate.tol

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
rcmp = function(n, lambda, nu, control = NULL)
{
	prep = prep.zicmp(n, lambda, nu)

	if (is.null(control)) {
		control = getOption("COMPoissonReg.control", default = get.control())
	}
	stopifnot("COMPoissonReg.control" %in% class(control))
	ymax = control$ymax
	hybrid.tol = control$hybrid.tol
	truncate.tol = control$truncate.tol

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
pcmp = function(x, lambda, nu, control = NULL)
{
	n = length(x)
	prep = prep.zicmp(n, lambda, nu)

	if (is.null(control)) {
		control = getOption("COMPoissonReg.control", default = get.control())
	}
	stopifnot("COMPoissonReg.control" %in% class(control))
	ymax = control$ymax
	hybrid.tol = control$hybrid.tol
	truncate.tol = control$truncate.tol

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
qcmp = function(q, lambda, nu, log.p = FALSE, control = NULL)
{
	n = length(q)
	prep = prep.zicmp(n, lambda, nu)

	if (is.null(control)) {
		control = getOption("COMPoissonReg.control", default = get.control())
	}
	stopifnot("COMPoissonReg.control" %in% class(control))
	ymax = control$ymax
	hybrid.tol = control$hybrid.tol
	truncate.tol = control$truncate.tol

	if (log.p) { lq = q } else { lq = log(q) }

	# We do this check here instead of in the C++ code, because we only want
	# a warning if quantiles were explicitly requested. (Not draws via rcmp,
	# for example).
	idx_warn = which(lq > log1p(-truncate.tol))
	if (length(idx_warn) > 0) {
		msg = sprintf(paste(
			"At least one requested quantile was very close to 1. In",
			"particular, %d of the given probabilities were greater than",
			"1 - truncate_tol = exp(%g), where truncate_tol = %g.",
			"Associated results may be affected by truncation. Consider",
			"adjusting the controls ymax and truncate.tol or reducing logq."),
			length(idx_warn), log1p(-truncate.tol), truncate.tol)
		warning(msg)
	}

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
ecmp = function(lambda, nu, control = NULL)
{
	# If lambda and nu are vectors, assume they are not repeats of the same
	# value and do not attempt to save time if this is the case.

	n = max(length(lambda), length(nu))
	prep = prep.zicmp(n, lambda, nu)

	if (is.null(control)) {
		control = getOption("COMPoissonReg.control", default = get.control())
	}
	stopifnot("COMPoissonReg.control" %in% class(control))
	ymax = control$ymax
	hybrid.tol = control$hybrid.tol
	truncate.tol = control$truncate.tol

	# If the following is true, use asymptotic approximation. Otherwise use
	# truncated sum.
	is.approx.valid = -1/prep$nu * log(prep$lambda) < log(hybrid.tol)

	out = numeric(n)
	for (i in 1:n) {
		if (is.approx.valid[i]) {
			# Use 1st derivative of the log-normalizing constant to compute the
			# expected value.
			out[i] = prep$lambda[i] * grad(ncmp, prep$lambda[i],
				nu = prep$nu[i], log = TRUE, control = control)
		} else {
			# Compute the expected value by a simple truncated sum
			x.seq = seq(0, tcmp(prep$lambda[i], prep$nu[i], control = control))
			out[i] = sum(x.seq * dcmp(x.seq, prep$lambda[i], prep$nu[i], control = control))
		}
	}

	return(out)
}

#' @name CMP Distribution
#' @export
vcmp = function(lambda, nu, control = NULL)
{
	# If lambda and nu are vectors, assume they are not repeats of the same
	# value and do not attempt to save time if this is the case.
	n = max(length(lambda), length(nu))
	prep = prep.zicmp(n, lambda, nu)

	if (is.null(control)) {
		control = getOption("COMPoissonReg.control", default = get.control())
	}
	stopifnot("COMPoissonReg.control" %in% class(control))
	ymax = control$ymax
	hybrid.tol = control$hybrid.tol
	truncate.tol = control$truncate.tol

	# If the following is true, use asymptotic approximation. Otherwise use
	# truncated sum.
	is.approx.valid = -1/prep$nu * log(prep$lambda) < log(hybrid.tol)

	# Compute expected value
	ev = ecmp(lambda, nu, control = control)

	out = numeric(n)
	for (i in 1:n) {
		if (is.approx.valid[i]) {
			# Use 2nd derivative of the log-normalizing constant to compute the
			# variance. The expression for dd is equal to Var(X) - E(X).
			dd = prep$lambda[i]^2 * hessian(ncmp, prep$lambda[i],
				nu = prep$nu[i], log = TRUE, control = control)
			out[i] = dd + ev[i]
		} else {
			# Compute the expected value by a simple truncated sum
			x.seq = seq(0, tcmp(prep$lambda[i], prep$nu[i], control = control))
			f.seq = dcmp(x.seq, prep$lambda[i], prep$nu[i], control = control)
			out[i] = sum(x.seq^2 * f.seq) - ev[i]^2
		}
	}

	return(out)
}

#' @name CMP Distribution
#' @export
ncmp = function(lambda, nu, log = FALSE, control = NULL)
{
	# If lambda and nu are vectors, assume they are not repeats of the same
	# value and do not attempt to save time if this is the case.
	n = max(length(lambda), length(nu))
	prep = prep.zicmp(n, lambda, nu)

	if (is.null(control)) {
		control = getOption("COMPoissonReg.control", default = get.control())
	}
	stopifnot("COMPoissonReg.control" %in% class(control))
	ymax = control$ymax
	hybrid.tol = control$hybrid.tol
	truncate.tol = control$truncate.tol

	out = z_hybrid(prep$lambda, prep$nu, take_log = log,
		hybrid_tol = hybrid.tol, truncate_tol = truncate.tol, ymax = ymax)

	return(out)
}

#' @name CMP Distribution
#' @export
tcmp = function(lambda, nu, control = NULL)
{
	if (is.null(control)) {
		control = getOption("COMPoissonReg.control", default = get.control())
	}
	stopifnot("COMPoissonReg.control" %in% class(control))
	ymax = control$ymax
	truncate.tol = control$truncate.tol

	n = max(length(lambda), length(nu))
	prep = prep.zicmp(n, lambda, nu)
	y_trunc(prep$lambda, prep$nu, tol = truncate.tol, ymax = ymax)

}
