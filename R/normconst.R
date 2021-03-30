#' CMP Normalizing Constant
#' 
#' @param lambda rate parameter.
#' @param nu dispersion parameter.
#' @param log If \code{TRUE}, return value on the log-scale.
#' @param method \code{trunc}, \code{approx}, or \code{hybrid}. See details.
#'
#' @details 
#' Compute the normalizing constant
#' \deqn{
#' Z(\lambda, \nu) = \sum_{x=0}^\infty \lambda^x / (x!)^\nu
#' }
#' 
#' Available methods are:
#' \describe{
#' \item{trunc}{Compute a truncated version of the series \eqn{Z(\lambda, \nu)}.}
#' \item{approx}{Use the approximation
#' \deqn{
#' Z(\lambda, \nu) = \frac{ \exp(\nu \lambda^{1/\nu}) }{ \lambda^{(\nu-1)/2\nu} (
#' 2\pi)^{(\nu-1)/2} \nu^{1/2} }
#' }
#' }
#' \item{hybrid}{Use the \code{approx} method if \eqn{\lambda^{1/\nu}} is small;
#' otherwise, use the \code{trunc} method. The \code{hybrid} method is used
#' when computing the density.}
#' }
#' 
#' @name CMP Normalizing Constant
#' @export
ncmp = function(lambda, nu, log = FALSE, method = "hybrid")
{
	ymax = getOption("COMPoissonReg.ymax", default = 1e6)
	n = max(length(lambda), length(nu))
	prep = prep.zicmp(n, lambda, nu)

	# These tolerances will be made into options in a later version of the
	# code. For now, too many places need to be changed from using default
	# values.
	hybrid_tol = 1e-2
	truncate_tol = 1e-6

	if (method == "trunc") {
		out = z_trunc(prep$lambda, prep$nu, tol = truncate_tol, take_log = log,
			ymax = ymax)
	} else if (method == "approx") {
		out = z_approx(prep$lambda, prep$nu, take_log = log)
	} else if (method == "hybrid") {
		out = z_hybrid(prep$lambda, prep$nu, take_log = log,
			tol1 = hybrid_tol, tol2 = truncate_tol)
	} else {
		msg = sprintf("Method not recognized: %s", method)
		stop(msg)
	}

	return(out)
}
