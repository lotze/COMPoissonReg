.onAttach = function(libname, pkgname) {
	options(COMPoissonReg.control = get.control())
	# options(COMPoissonReg.ymax = 1e6)
	# options(COMPoissonReg.optim.method = 'L-BFGS-B')
	# options(COMPoissonReg.optim.control = list(maxit = 150))
	# options(COMPoissonReg.grad.eps = 1e-5)
	# options(COMPoissonReg.hess.eps = 1e-2)
	# options(COMPoissonReg.hybrid.tol = 1e-2)
	# options(COMPoissonReg.truncate.tol = 1e-6)
}

# TBD: Needs documentation
#' @export
get.fixed = function(beta = integer(0L), gamma = integer(0L), zeta = integer(0L))
{
	stopifnot(is.integer(beta))
	stopifnot(is.integer(gamma))
	stopifnot(is.integer(zeta))
	out = list(
		beta = sort(unique(beta)),
		gamma = sort(unique(gamma)),
		zeta = sort(unique(zeta))
	)
	class(out) = "glm.cmp.fixed"
	return(out)
}

# TBD: Needs documentation
#' @export
get.init = function(beta = NULL, gamma = NULL, zeta = NULL)
{
	out = list(beta = beta, gamma = gamma, zeta = zeta)
	class(out) = "glm.cmp.init"
	return(out)
}

# TBD: Needs documentation
#' @export
get.offset = function(x = 0, s = 0, w = 0)
{
	out = list(x = x, s = s, w = s)
	class(out) = "glm.cmp.offset"
	return(out)
}

# TBD: Needs documentation
#' @param z.method a string: \code{hybrid}, \code{approx}, or \code{trunc}. TBD: explain
#' @export
get.control = function(ymax = 1e6, z.method = "hybrid", optim.method = 'L-BFGS-B',
	optim.control = list(maxit = 150), grad.eps = 1e-5, hess.eps = 1e-2,
	hybrid.tol = 1e-2, truncate.tol = 1e-6)
{
	out = list(ymax = ymax, z.method = z.method, optim.method = optim.method,
		optim.control = optim.control, grad.eps = grad.eps, hess.eps = hess.eps,
		hybrid.tol = hybrid.tol, truncate.tol = truncate.tol)
	class(out) = "COMPoissonReg.control"
	return(out)
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

format.difftime = function(x) {
	s = as.numeric(x, units = "secs")
	dd = floor(s / (60^2 * 24))
	dd.resid = s / (60^2 * 24) - dd
	hh = floor(24*dd.resid)
	hh.resid = 24*dd.resid - floor(24*dd.resid)
	mm = floor(60*hh.resid)
	mm.resid = 60*hh.resid - floor(60*hh.resid)
	ss = floor(60*mm.resid)

	if (dd > 0) {
		fmt = sprintf("%02dd:%02dh:%02dm:%02ds", dd, hh, mm, ss)
	} else if (hh > 0) {
		fmt = sprintf("%02dh:%02dm:%02ds", hh, mm, ss)
	} else if (mm > 0) {
		fmt = sprintf("%02dm:%02ds", mm, ss)
	} else {
		fmt = sprintf("%d sec", ss)
	}

	return(fmt)
}

printf = function(msg, ...) {
	cat(sprintf(msg, ...))
}

logger = function(msg, ...)
{
	sys.time = as.character(Sys.time())
	cat(sys.time, "-", sprintf(msg, ...))
}

grad.fwd = function(f, x, h = 1e-5, ...) {
	k = length(x)
	eye = diag(1, k)
	res = numeric(k)
	fx = f(x, ...)
	for (j in 1:k) {
		res[j] = ( f(x + h * eye[,j], ...) - fx ) / h
	}
	return(res)
}

hess.fwd = function(f, x, h = 1e-5, ...) {
	k = length(x)
	eye = diag(1, k)
	H = matrix(NA, k, k)

	fx = f(x, ...)
	fx.eps = numeric(k)
	for (j in 1:k) {
		fx.eps[j] = f(x + h * eye[,j], ...)
	}

	for (j in 1:k) {
		for (l in 1:k) {
			num = f(x + h * eye[,j] + h * eye[,l], ...) -
				fx.eps[l] - fx.eps[j] + fx
			H[j,l] = num / h^2
		}
	}
	(H + t(H)) / 2
}

is.zero.matrix = function(X, eps = 1e-12)
{
	all(abs(X) < eps)
}

is.intercept.only = function(X, eps = 1e-12)
{
	n = length(X)
	all(dim(X) == c(n,1)) & is.zero.matrix(X-1, eps = eps)
}
