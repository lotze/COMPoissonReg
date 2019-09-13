fit_zicmp_reg = function(y, X, S, W, beta_init, gamma_init, zeta_init, off_x, off_s, off_w)
{
	start = Sys.time()
	u = as.integer(y == 0)
	d1 = ncol(X)
	d2 = ncol(S)
	d3 = ncol(W)
	n = length(y)
	qq = d1 + d2 + d3

	if (is.null(beta_init)) { beta_init = rep(0, d1) }
	if (is.null(gamma_init)) { gamma_init = rep(0, d2) }
	if (is.null(zeta_init)) { zeta_init = rep(0, d3) }
	stopifnot(d1 == length(beta_init))
	stopifnot(d2 == length(gamma_init))
	stopifnot(d3 == length(zeta_init))

	optim_method = getOption("COMPoissonReg.optim.method")
	optim_control = getOption("COMPoissonReg.optim.control")

	tx = function(par) {
		list(
			beta = par[1:d1],
			gamma = par[1:d2 + d1],
			zeta = par[1:d3 + d1 + d2]
		)
	}

	loglik = function(par) {
		theta = tx(par)
		out = fitted_zicmp_internal(X, S, W, theta$beta, theta$gamma, theta$zeta, off_x, off_s, off_w)
		lambda = out$lambda
		nu = out$nu
		p = out$p
		logz = z_hybrid(lambda, nu, take_log = TRUE)
		t(u) %*% log(p*exp(logz) + (1-p)) + t(1 - u) %*% (log(1-p) +
			y*log(lambda) - nu*lgamma(y+1)) - sum(logz)
	}

	optim_control$fnscale = -1
	par_init = c(beta_init, gamma_init, zeta_init)
	res = optim(par_init, loglik, method = optim_method,
		control = optim_control, hessian = TRUE)

	theta_hat = list(
		beta = res$par[1:d1],
		gamma = res$par[1:d2 + d1],
		zeta = res$par[1:d3 + d1 + d2]
	)
	names(theta_hat$beta) = sprintf("X:%s", colnames(X))
	names(theta_hat$gamma) = sprintf("S:%s", colnames(S))
	names(theta_hat$zeta) = sprintf("W:%s", colnames(W))

	H = res$hessian
	colnames(H) = rownames(H) = c(
		sprintf("X:%s", colnames(X)),
		sprintf("S:%s", colnames(S)),
		sprintf("W:%s", colnames(W)))

	loglik = res$value
	elapsed_sec = as.numeric(Sys.time() - start, type = "sec")

	res = list(theta_hat = theta_hat, H = H, opt_res = res,
		elapsed_sec = elapsed_sec, loglik = loglik, n = n)
	return(res)
}

fit_cmp_reg = function(y, X, S, beta_init, gamma_init, off_x, off_s)
{
	start = Sys.time()
	u = as.integer(y == 0)
	d1 = ncol(X)
	d2 = ncol(S)
	n = length(y)
	qq = d1 + d2

	if (is.null(beta_init)) { beta_init = rep(0, d1) }
	if (is.null(gamma_init)) { gamma_init = rep(0, d2) }
	stopifnot(d1 == length(beta_init))
	stopifnot(d2 == length(gamma_init))

	optim_method = getOption("COMPoissonReg.optim.method")
	optim_control = getOption("COMPoissonReg.optim.control")

	tx = function(par) {
		list(
			beta = par[1:d1],
			gamma = par[1:d2 + d1]
		)
	}

	loglik = function(par) {
		theta = tx(par)
		out = fitted_cmp_internal(X, S, theta$beta, theta$gamma, off_x, off_s)
		lambda = out$lambda
		nu = out$nu

		logz = z_hybrid(lambda, nu, take_log = TRUE)
		sum(y*log(lambda) - nu*lgamma(y+1) - logz)
	}

	optim_control$fnscale = -1
	par_init = c(beta_init, gamma_init)

	res = optim(par_init, loglik, method = optim_method,
		control = optim_control, hessian = TRUE)

	theta_hat = list(
		beta = res$par[1:d1],
		gamma = res$par[1:d2 + d1]
	)
	names(theta_hat$beta) = sprintf("X:%s", colnames(X))
	names(theta_hat$gamma) = sprintf("S:%s", colnames(S))

	H = res$hessian
	colnames(H) = rownames(H) = c(
		sprintf("X:%s", colnames(X)),
		sprintf("S:%s", colnames(S))
	)

	loglik = res$value
	elapsed_sec = as.numeric(Sys.time() - start, type = "sec")

	res = list(theta_hat = theta_hat, H = H, opt_res = res,
		elapsed_sec = elapsed_sec, loglik = loglik, n = n)
	return(res)
}

fit_zip_reg = function(y, X, W, beta_init, zeta_init, off_x, off_w)
{
	start = Sys.time()
	u = as.integer(y == 0)
	d1 = ncol(X)
	d3 = ncol(W)
	n = length(y)
	qq = d1 + d3

	if (is.null(beta_init)) { beta_init = rep(0, d1) }
	if (is.null(zeta_init)) { zeta_init = rep(0, d3) }
	stopifnot(d1 == length(beta_init))
	stopifnot(d3 == length(zeta_init))

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
		out = fitted_zip_internal(X, W, theta$beta, theta$zeta, off_x, off_w)
		lambda = out$lambda
		p = out$p

		sum(u*log(p + (1-p)*exp(-lambda)) + (1-u)*log(1-p) +
			(1-u)*(y*log(lambda) - lambda - lgamma(y+1)))
	}

	optim_control$fnscale = -1
	par_init = c(beta_init, zeta_init)
	res = optim(par_init, loglik, method = optim_method,
		control = optim_control, hessian = TRUE)

	theta_hat = list(
		beta = res$par[1:d1],
		zeta = res$par[1:d3 + d1]
	)

	H = res$hessian
	colnames(H) = rownames(H) = c(
		sprintf("X:%s", colnames(X)),
		sprintf("W:%s", colnames(W)))

	loglik = res$value
	elapsed_sec = as.numeric(Sys.time() - start, type = "sec")

	res = list(theta_hat = theta_hat, H = H, opt_res = res,
		elapsed_sec = elapsed_sec, loglik = loglik, n = n)
	return(res)
}
