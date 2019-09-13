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
#' @param log_p logical; if TRUE, probabilities \code{p} are given as \eqn{\log(p)}.
#' 
#' @return 
#' \itemize{
#' \item \code{dcmp} gives the density,
#' \item \code{pcmp} gives the cumulative probability,
#' \item \code{qcmp} gives the quantile function, and
#' \item \code{rcmp} generates random values.
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
	prep = prep_zicmp(length(x), lambda, nu)
	dcmp_cpp(x, prep$lambda, prep$nu, take_log = log)
}

#' @name CMP Distribution
#' @export
rcmp = function(n, lambda, nu)
{
	prep = prep_zicmp(n, lambda, nu)
	ymax = getOption("COMPoissonReg.ymax")
	rcmp_cpp(n, prep$lambda, prep$nu, ymax = ymax)
}

#' @name CMP Distribution
#' @export
pcmp = function(x, lambda, nu)
{
	prep = prep_zicmp(length(x), lambda, nu)
	pcmp_cpp(x, prep$lambda, prep$nu)
}

#' @name CMP Distribution
#' @export
qcmp = function(q, lambda, nu, log_p = FALSE)
{
	prep = prep_zicmp(length(q), lambda, nu)
	if (log_p) {
		log_q = q
	} else {
		log_q = log(q)
	}
	ymax = getOption("COMPoissonReg.ymax")
	qcmp_cpp(log_q, prep$lambda, prep$nu, ymax = ymax)
}

cmp_expected_value = function(lambda, nu)
{
	n = max(length(lambda), length(nu))
	prep = prep_zicmp(n, lambda, nu)

	res = numeric(n)
	for (i in 1:n) {
		res[i] = prep$lambda[i] * grad_fwd(z_hybrid, prep$lambda[i], nu = prep$nu[i], take_log = TRUE)
	}
	return(res)
}
