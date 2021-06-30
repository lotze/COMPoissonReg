fit.zicmp.reg = function(y, X, S, W, beta.init, gamma.init, zeta.init,
	off.x, off.s, off.w, fixed.beta, fixed.gamma, fixed.zeta)
{
	start = Sys.time()
	u = as.integer(y == 0)
	d1 = ncol(X)
	d2 = ncol(S)
	d3 = ncol(W)
	n = length(y)
	qq = d1 + d2 + d3

	if (is.null(beta.init)) { beta.init = numeric(d1) }
	if (is.null(gamma.init)) { gamma.init = numeric(d2) }
	if (is.null(zeta.init)) { zeta.init = numeric(d3) }
	stopifnot(d1 == length(beta.init))
	stopifnot(d2 == length(gamma.init))
	stopifnot(d3 == length(zeta.init))

	optim.method = getOption("COMPoissonReg.optim.method")
	optim.control = getOption("COMPoissonReg.optim.control")
	hybrid.tol = getOption("COMPoissonReg.hybrid.tol")
	truncate.tol = getOption("COMPoissonReg.truncate.tol")
	ymax = getOption("COMPoissonReg.ymax")

	# "par" represents the vector of optimization variables; it contains only
	# coefficients which have not been fixed.
	#
	# "theta" is a list whose three elements represent the three vectors of
	# coefficients; these include fixed variables in addition to optimization
	# variables. Fixed variables are assumed to be set to their initial values.

	fixed.beta = sort(unique(fixed.beta))
	fixed.gamma = sort(unique(fixed.gamma))
	fixed.zeta = sort(unique(fixed.zeta))
	unfixed.beta = sort(setdiff(seq_len(d1), fixed.beta))
	unfixed.gamma = sort(setdiff(seq_len(d2), fixed.gamma))
	unfixed.zeta = sort(setdiff(seq_len(d3), fixed.zeta))
	idx.par1 = seq_along(unfixed.beta)
	idx.par2 = seq_along(unfixed.gamma) + length(unfixed.beta)
	idx.par3 = seq_along(unfixed.zeta) + length(unfixed.beta) +  length(unfixed.gamma)

	# The following functions transform back and forth between the "par" and
	# "theta" representations.
	
	par2theta = function(par) {
		beta = rep(NA, d1)
		beta[fixed.beta] = beta.init[fixed.beta]
		beta[unfixed.beta] = par[idx.par1]

		gamma = rep(NA, d2)
		gamma[fixed.gamma] = gamma.init[fixed.gamma]
		gamma[unfixed.gamma] = par[idx.par2]

		zeta = rep(NA, d3)
		zeta[fixed.zeta] = zeta.init[fixed.zeta]
		zeta[unfixed.zeta] = par[idx.par3]

		list(beta = beta, gamma = gamma, zeta = zeta)
	}

	theta2par = function(theta) {
		c(theta$beta[unfixed.beta], theta$gamma[unfixed.gamma], theta$zeta[unfixed.zeta])
	}

	loglik = function(par) {
		theta = par2theta(par)
		out = fitted.zicmp.internal(X, S, W, theta$beta, theta$gamma,
			theta$zeta, off.x, off.s, off.w)
		loglik_zicmp(y, out$lambda, out$nu, out$p, hybrid.tol, truncate.tol, ymax)
	}

	if (!is.null(optim.control$fnscale)) {
		warning("COMPoissonReg.optim.control$fnscale disregarded and taken as -1")
	}
	optim.control$fnscale = -1
	par.init = theta2par(list(beta = beta.init, gamma = gamma.init, zeta = zeta.init))
	res = optim(par.init, loglik, method = optim.method,
		control = optim.control, hessian = TRUE)

	# theta includes fixed values as well as optimization results.
	
	theta.hat = par2theta(res$par)
	names(theta.hat$beta) = sprintf("X:%s", colnames(X))
	names(theta.hat$gamma) = sprintf("S:%s", colnames(S))
	names(theta.hat$zeta) = sprintf("W:%s", colnames(W))

	# The Hessian from the optimizer only has entries corresponding to the
	# optimization variables. To apply labels from the design matrices, we
	# must pick out the non-fixed variables from each.
	
	H = res$hessian
	colnames(H) = rownames(H) = c(
		sprintf("X:%s", colnames(X)[unfixed.beta]),
		sprintf("S:%s", colnames(S)[unfixed.gamma]),
		sprintf("W:%s", colnames(W)[unfixed.zeta]))

	loglik = res$value
	elapsed.sec = as.numeric(Sys.time() - start, units = "secs")

	res = list(theta.hat = theta.hat, H = H, opt.res = res,
		elapsed.sec = elapsed.sec, loglik = loglik, n = n,
		optim.method = optim.method, optim.control = optim.control,
		fixed.beta = fixed.beta, fixed.gamma = fixed.gamma,
		fixed.zeta = fixed.zeta, unfixed.beta = unfixed.beta,
		unfixed.gamma = unfixed.gamma, unfixed.zeta = unfixed.zeta)
	return(res)
}

# This is just provided as a convenience to call fit.zicmp.reg with some dummy
# values.
fit.cmp.reg = function(y, X, S, beta.init, gamma.init, off.x, off.s, fixed.beta,
	fixed.gamma)
{
	n = length(y)
	W = matrix(NA, n, 0)
	zeta.init = numeric(0)
	off.w = numeric(n)
	fixed.zeta = integer(0)

	fit.zicmp.reg(y, X, S, W, beta.init, gamma.init, zeta.init,
		off.x, off.s, off.w, fixed.beta, fixed.gamma, fixed.zeta)
}
