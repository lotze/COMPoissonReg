.onLoad <- function(libname, pkgname){
	options(COMPoissonReg.optim.method = 'L-BFGS-B')
	options(COMPoissonReg.optim.control = list(maxit = 150))
}

logger <- function(msg, ...)
{
	sys.time <- as.character(Sys.time())
	cat(sys.time, "-", sprintf(msg, ...))
}
