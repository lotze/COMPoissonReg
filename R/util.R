.onAttach = function(libname, pkgname){
	options(COMPoissonReg.optim.method = 'L-BFGS-B')
	options(COMPoissonReg.optim.control = list(maxit = 150))
	options(COMPoissonReg.grad.eps = 1e-5)
	options(COMPoissonReg.hess.eps = 1e-2)
	options(COMPoissonReg.ymax = 1e6)
}

format_difftime = function(x) {
	s = as.numeric(x, units = "secs")
	dd = floor(s / (60^2 * 24))
	dd_resid = s / (60^2 * 24) - dd
	hh = floor(24*dd_resid)
	hh_resid = 24*dd_resid - floor(24*dd_resid)
	mm = floor(60*hh_resid)
	mm_resid = 60*hh_resid - floor(60*hh_resid)
	ss = floor(60*mm_resid)

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
