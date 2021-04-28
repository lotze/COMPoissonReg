fit.zicmp.reg = function(y, X, S, W, beta.init, gamma.init, zeta.init, off.x, off.s, off.w)
{
	start = Sys.time()
	u = as.integer(y == 0)
	d1 = ncol(X)
	d2 = ncol(S)
	d3 = ncol(W)
	n = length(y)
	qq = d1 + d2 + d3

	if (is.null(beta.init)) { beta.init = rep(0, d1) }
	if (is.null(gamma.init)) { gamma.init = rep(0, d2) }
	if (is.null(zeta.init)) { zeta.init = rep(0, d3) }
	stopifnot(d1 == length(beta.init))
	stopifnot(d2 == length(gamma.init))
	stopifnot(d3 == length(zeta.init))

	optim_method = getOption("COMPoissonReg.optim.method")
	optim_control = getOption("COMPoissonReg.optim.control")
	hybrid_tol = getOption("COMPoissonReg.hybrid.tol")
	truncate_tol = getOption("COMPoissonReg.truncate.tol")
	ymax = getOption("COMPoissonReg.ymax")

	tx = function(par) {
		list(
			beta = par[1:d1],
			gamma = par[1:d2 + d1],
			zeta = par[1:d3 + d1 + d2]
		)
	}

	loglik = function(par) {
		theta = tx(par)
		out = fitted_zicmp_internal(X, S, W, theta$beta, theta$gamma,
			theta$zeta, off.x, off.s, off.w)
		loglik_zicmp(y, out$lambda, out$nu, out$p, hybrid_tol, truncate_tol,
			ymax);
	}

	optim_control$fnscale = -1
	par.init = c(beta.init, gamma.init, zeta.init)
	res = optim(par.init, loglik, method = optim_method,
		control = optim_control, hessian = TRUE)

	theta.hat = list(
		beta = res$par[1:d1],
		gamma = res$par[1:d2 + d1],
		zeta = res$par[1:d3 + d1 + d2]
	)
	names(theta.hat$beta) = sprintf("X:%s", colnames(X))
	names(theta.hat$gamma) = sprintf("S:%s", colnames(S))
	names(theta.hat$zeta) = sprintf("W:%s", colnames(W))

	H = res$hessian
	colnames(H) = rownames(H) = c(
		sprintf("X:%s", colnames(X)),
		sprintf("S:%s", colnames(S)),
		sprintf("W:%s", colnames(W)))

	loglik = res$value
	elapsed.sec = as.numeric(Sys.time() - start, units = "secs")

	res = list(theta.hat = theta.hat, H = H, opt.res = res,
		elapsed.sec = elapsed.sec, loglik = loglik, n = n)
	return(res)
}

fit.cmp.reg = function(y, X, S, beta.init, gamma.init, off.x, off.s)
{
	start = Sys.time()
	u = as.integer(y == 0)
	d1 = ncol(X)
	d2 = ncol(S)
	n = length(y)
	qq = d1 + d2

	if (is.null(beta.init)) { beta.init = rep(0, d1) }
	if (is.null(gamma.init)) { gamma.init = rep(0, d2) }
	stopifnot(d1 == length(beta.init))
	stopifnot(d2 == length(gamma.init))

	optim_method = getOption("COMPoissonReg.optim.method")
	optim_control = getOption("COMPoissonReg.optim.control")
	hybrid_tol = getOption("COMPoissonReg.hybrid.tol")
	truncate_tol = getOption("COMPoissonReg.truncate.tol")
	ymax = getOption("COMPoissonReg.ymax")

	tx = function(par) {
		list(
			beta = par[1:d1],
			gamma = par[1:d2 + d1]
		)
	}

	loglik = function(par) {
		theta = tx(par)
		out = fitted_cmp_internal(X, S, theta$beta, theta$gamma, off.x, off.s)
		loglik_cmp(y, out$lambda, out$nu, hybrid_tol, truncate_tol, ymax)
	}

	optim_control$fnscale = -1
	par.init = c(beta.init, gamma.init)

	res = optim(par.init, loglik, method = optim_method,
		control = optim_control, hessian = TRUE)

	theta.hat = list(
		beta = res$par[1:d1],
		gamma = res$par[1:d2 + d1]
	)
	names(theta.hat$beta) = sprintf("X:%s", colnames(X))
	names(theta.hat$gamma) = sprintf("S:%s", colnames(S))

	H = res$hessian
	colnames(H) = rownames(H) = c(
		sprintf("X:%s", colnames(X)),
		sprintf("S:%s", colnames(S))
	)

	loglik = res$value
	elapsed.sec = as.numeric(Sys.time() - start, units = "secs")

	res = list(theta.hat = theta.hat, H = H, opt.res = res,
		elapsed.sec = elapsed.sec, loglik = loglik, n = n)
	return(res)
}

fit.zip.reg = function(y, X, W, beta.init, zeta.init, off.x, off.w)
{
	start = Sys.time()
	u = as.integer(y == 0)
	d1 = ncol(X)
	d3 = ncol(W)
	n = length(y)
	qq = d1 + d3

	if (is.null(beta.init)) { beta.init = rep(0, d1) }
	if (is.null(zeta.init)) { zeta.init = rep(0, d3) }
	stopifnot(d1 == length(beta.init))
	stopifnot(d3 == length(zeta.init))

	optim_method = getOption("COMPoissonReg.optim.method")
	optim_control = getOption("COMPoissonReg.optim.control")

	tx = function(par) {
		list(
			beta = par[1:d1],
			zeta = par[1:d3 + d1]
		)
	}

	loglik = function(par) {
		theta = tx(par)
		out = fitted_zip_internal(X, W, theta$beta, theta$zeta, off.x, off.w)
		lambda = out$lambda
		p = out$p

		sum(u*log(p + (1-p)*exp(-lambda)) + (1-u)*log(1-p) +
			(1-u)*(y*log(lambda) - lambda - lgamma(y+1)))
	}

	optim_control$fnscale = -1
	par.init = c(beta.init, zeta.init)
	res = optim(par.init, loglik, method = optim_method,
		control = optim_control, hessian = TRUE)

	theta.hat = list(
		beta = res$par[1:d1],
		zeta = res$par[1:d3 + d1]
	)

	H = res$hessian
	colnames(H) = rownames(H) = c(
		sprintf("X:%s", colnames(X)),
		sprintf("W:%s", colnames(W)))

	loglik = res$value
	elapsed.sec = as.numeric(Sys.time() - start, units = "secs")

	res = list(theta.hat = theta.hat, H = H, opt.res = res,
		elapsed.sec = elapsed.sec, loglik = loglik, n = n)
	return(res)
}

