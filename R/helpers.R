.onAttach = function(libname, pkgname) {
	options(COMPoissonReg.control = get.control())
}

#' Construct an object that specifies which indices of coefficients should
#' remain fixed in maximum likelihood computation.
#'
#' @param beta Vector of indices of \code{beta} to keep fixed.
#' @param gamma Vector of indices of \code{gamma} to keep fixed.
#' @param zeta Vector of indices of \code{zeta} to keep fixed.
#'
#' @details
#' Arguments are expected to be vectors of integers. These are interpreted as
#' the indices to keep fixed during optimization. For example,
#' \code{beta = c(1L, 1L, 2L)} indicates that the first and second elements of
#' \code{beta} should remain fixed. Note that duplicate indices are ignored.
#' The default value is the empty vector \code{integer(0)}, which requests that
#' no elements of the given coefficient vector should be fixed. 
#'
#' @return List of vectors indicating fixed indices.
#' @export
get.fixed = function(beta = integer(0), gamma = integer(0), zeta = integer(0))
{
	stopifnot(is.integer(beta))
	stopifnot(is.integer(gamma))
	stopifnot(is.integer(zeta))
	out = list(
		beta = sort(unique(beta)),
		gamma = sort(unique(gamma)),
		zeta = sort(unique(zeta))
	)
	class(out) = "COMPoissonReg.fixed"
	return(out)
}

#' Construct initial values for coefficients with zeros.
#'
#' @param d1 Dimension of \code{beta}.
#' @param d2 Dimension of \code{gamma}.
#' @param d3 Dimension of \code{zeta}.
#'
#' @return List of initial value terms containing all zeros.
#' @export
get.init.zero = function(d1 = 0, d2 = 0, d3 = 0)
{
	get.init(beta = numeric(d1), gamma = numeric(d2), zeta = numeric(d3))
}

#' Construct initial values for coefficients.
#'
#' @param beta Vector for \code{beta}.
#' @param gamma Vector for \code{gamma}.
#' @param zeta Vector for \code{zeta}.
#'
#' @details
#' The default value \code{NULL} is interpreted as an empty vector, so that the
#' given component is absent from the model.
#'
#' @return List of initial value terms.
#' @export
get.init = function(beta = NULL, gamma = NULL, zeta = NULL)
{
	if (is.null(beta)) { beta = numeric(0) }
	if (is.null(gamma)) { gamma = numeric(0) }
	if (is.null(zeta)) { zeta = numeric(0) }
	out = list(beta = beta, gamma = gamma, zeta = zeta)
	class(out) = "COMPoissonReg.init"
	return(out)
}

#' Construct zero values for offsets.
#'
#' @param n Number of observations.
#'
#' @return List of offset terms containing all zeros.
#' @export
get.offset.zero = function(n)
{
	get.offset(x = numeric(n), s = numeric(n), w = numeric(n))
}

#' Construct values for offsets.
#'
#' @param x Vector of offsets to go with \code{X} matrix.
#' @param s Vector of offsets to go with \code{S} matrix.
#' @param w Vector of offsets to go with \code{W} matrix.
#'
#' @details
#' The default value \code{NULL} is interpreted as a vector of zeros. At least
#' one component must be non-NULL so that the dimension can be determined.
#'
#' @return List of offset terms.
#' @export
get.offset = function(x = NULL, s = NULL, w = NULL)
{
	lengths = c(length(x), length(s), length(w))

	# At least one offset should be non-null. Set n to be the first entry
	stopifnot(!is.null(lengths))
	n = max(lengths)

	if (is.null(x)) { x = numeric(n) }
	if (is.null(s)) { s = numeric(n) }
	if (is.null(w)) { w = numeric(n) }
	stopifnot(n == length(x))
	stopifnot(n == length(s))
	stopifnot(n == length(w))

	out = list(x = x, s = s, w = s)
	class(out) = "COMPoissonReg.offset"
	return(out)
}

#' Construct a control object to pass additional arguments to a number of
#' functions in the package.
#' 
#' @param ymax Truncate counts to maximum value of \code{y}.
#' @param z.method One of: \code{"hybrid"}, \code{"approx"}, or \code{"trunc"}.
#' See details.
#' @param optim.method Optimization method for maximum likelihood. See the
#' \code{method} argument in \link[stats]{optim}.
#' @param optim.control \code{control} argument for \link[stats]{optim}.
#' @param hybrid.tol Tolerance for \code{z.method = "hybrid"}. See details.
#' @param truncate.tol Tolerance for \code{z.method = "trunc"}. See details.
#' 
#' @details
#' The element \code{z.method} specifies how the CMP normalizing constant and
#' related quantities are computed.
#' \itemize{
#' \item \code{z.method = "approx"}: use an approximation which is easy to
#' compute but suitable when \eqn{\lambda^{-1/\nu}} is small.
#' 
#' \item \code{z.method = "trunc"}: the infinite series for the normalizing
#' constant is truncated based on \code{truncate.tol}. The procedure is
#' described in the package vignette.
#' 
#' \item \code{z.method = "hybrid"}: use the \code{approx} method when
#' \eqn{\lambda^{-1/\nu}} is smaller than \code{hybrid.tol}; otherwise use the
#' \code{trunc} method.
#' }
#' 
#' The element \code{ymax} protects against very long computations. Users
#' should beware when increasing this significantly beyond the default, as it
#' may result in a session which needs to be terminated.
#' 
#' @return List of controls.
#' @export
get.control = function(ymax = 1e6, z.method = "hybrid", optim.method = 'L-BFGS-B',
	optim.control = list(maxit = 150), hybrid.tol = 1e-2, truncate.tol = 1e-6)
{
	stopifnot(z.method %in% c("hybrid", "approx",  "trunc"))

	out = list(ymax = ymax, z.method = z.method, optim.method = optim.method,
		optim.control = optim.control, hybrid.tol = hybrid.tol,
		truncate.tol = truncate.tol)
	class(out) = "COMPoissonReg.control"
	return(out)
}

#' Construct model matrices and offsets for CMP/ZICMP regression
#'
#' @param X An \code{X} matrix to use with \code{beta}.
#' @param S An \code{S} matrix to use with \code{gamma}.
#' @param W A \code{W} matrix to use with \code{zeta}.
#' @param offset An offset object. See helper function \link{get.offset}.
#'
#' @return List of model matrix terms.
#' @export
get.modelmatrix = function(X = NULL, S = NULL, W = NULL, offset = NULL)
{
	nrows = c(nrow(X), nrow(S), nrow(W))

	# At least one design matrix should be non-null. Set n to be the first entry
	stopifnot(!is.null(nrows))
	n = head(nrows, 1)

	if (is.null(X)) { X = matrix(0, n, 0) }
	if (is.null(S)) { S = matrix(0, n, 0) }
	if (is.null(W)) { W = matrix(0, n, 0) }
	stopifnot(n == nrow(X))
	stopifnot(n == nrow(S))
	stopifnot(n == nrow(W))

	if (is.null(offset)) {
		offset = get.offset.zero(n)
	}
	stopifnot("COMPoissonReg.offset" %in% class(offset))

	out = list(X = X, S = S, W = W, offset = offset)
	class(out) = "COMPoissonReg.modelmatrix"
	return(out)
}

#' Prepare lambda, nu, and p in vector form for use with CMP/ZICMP distribution
#' functions.
#' 
#' @param n Number of observations
#' @param lambda The rate parameter: scalar or vector of length n
#' @param nu The dispersion parameter: scalar or vector of length n
#' @param p The zero-inflation parameter: scalar or vector of length n
#' 
#' @details 
#' Extend lambda, nu, and p vectors to be compatible lengths. If all are length
#' 1, do not extend them - this is a special case which is handled more
#' efficiently. Also make sure parameters are in the right space
#' 
#' @noRd
prep.zicmp = function(n, lambda, nu, p = 0)
{
	L = max(length(lambda), length(nu), length(p))

	stopifnot(all(lambda >= 0))
	stopifnot(all(nu >= 0))
	stopifnot(all(p >= 0 & p <= 1))

	if (n > 1 && L > 1) { stopifnot(n == L) }
	if (length(lambda) == 1 && L > 1) { lambda = rep(lambda, L) }
	if (length(nu) == 1 && L > 1) { nu = rep(nu, L) }
	if (length(p) == 1 && L > 1) { p = rep(p, L) }
	if (L > 1) { type = "indep" } else { type = "iid" }

	list(lambda = lambda, nu = nu, p = p, type = type)
}
