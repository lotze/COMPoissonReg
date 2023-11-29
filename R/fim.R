# These functions are currently not exported to users. Numerical computation
# of the FIM seems to be potentially unstable, whether computing by
# derivatives of the normalizing function, or explicitly programming the
# infinite sum expressions and truncating them. Instead, we currently
# provide variances via the Hessian from the optimizer.

# Compute the CMP information matrix using the exact expression, except
# use numerical derivatives of z functions. This avoids numerical issues
# with infinite sums.
fim.cmp = function(lambda, nu)
{
	f = function(theta) {
		z_hybrid(theta[1], theta[2], take_log = TRUE)
	}
	hess.eps = getOption("COMPoissonReg.hess.eps", default = 1e-2)
	H = hess_fwd(f, c(lambda, nu), h = hess.eps)
	H[1,1] = H[1,1] + ecmp(lambda, nu) / lambda^2
	return(H)
}

# Compute the CMP information matrix using Monte Carlo. This avoids
# issues with numerical second derivatives.
fim.cmp.mc = function(lambda, nu, reps)
{
	stopifnot(length(lambda) == 1)
	stopifnot(length(nu) == 1)
	grad.eps = getOption("COMPoissonReg.grad.eps", default = 1e-5)

	f = function(theta, data) {
		sum(dcmp(data, lambda = theta[1], nu = theta[2], log = TRUE))
	}

	FIM = matrix(0, 2, 2)
	x = rcmp(reps, lambda, nu)
	for (r in 1:reps) {
		S = grad_fwd(f, x = c(lambda, nu), h = grad.eps, data = x[r])
		FIM = FIM + S %*% t(S)
	}

	return(FIM / reps)
}

# Compute the ZICMP information matrix using the exact expression, except
# use numerical derivatives of z functions.
fim.zicmp = function(lambda, nu, p)
{
	stopifnot(length(lambda) == 1)
	stopifnot(length(nu) == 1)
	stopifnot(length(p) == 1)

	z = z_hybrid(lambda, nu)
	meany = ezicmp(lambda, nu, p)
	f = function(theta) {
		z_hybrid(theta[1], theta[2], take_log = TRUE)
	}

	FIM = matrix(NA, 3, 3)
	colnames(FIM) = c("lambda", "nu", "p")
	rownames(FIM) = c("lambda", "nu", "p")

	# The gradient and hessian of the z-like functions all have a closed
	# form, but require some work to compute reliably. We have some code to
	# compute the z-function itself somewhat reliably, so it's safer
	# to use numerical derivatives here. We use forward derivatives
	# to avoid issues at the boundary when lambda or nu is near 0.
	grad.eps = getOption("COMPoissonReg.grad.eps", default = 1e-5)
	hess.eps = getOption("COMPoissonReg.hess.eps", default = 1e-2)
	Dlogz = grad_fwd(f, c(lambda, nu), h = grad.eps)
	Hlogz = hess_fwd(f, c(lambda, nu), h = hess.eps)

	dlogzdlambda = Dlogz[1]
	dlogzdnu = Dlogz[2]
	d2logzdlambda2 = Hlogz[1,1]
	d2logzdnu2 = Hlogz[2,2]
	d2logzdlambdadnu = Hlogz[1,2]

	# FIM[lambda, lambda]
	FIM[1,1] = (1-p)*d2logzdlambda2 - p*(1-p)*dlogzdlambda^2 / (p*(z-1)+1) + meany / lambda^2

	# FIM[nu, nu]
	FIM[2,2] = (1-p)*d2logzdnu2 - p*(1-p)*dlogzdnu^2 / (p*(z-1)+1)

	# FIM[p, p]
	FIM[3,3] = (1/z) * (z-1)^2 / (p*(z-1) + 1) + (1 - 1/z) / (1-p)

	# FIM[lambda, nu]
	FIM[1,2] = (1-p)*d2logzdlambdadnu - p*(1-p)*dlogzdnu*dlogzdlambda / (p*(z-1)+1)
	FIM[2,1] = FIM[1,2]

	# FIM[lambda, p]
	FIM[1,3] = -1/(p*(z-1) + 1) * dlogzdlambda
	FIM[3,1] = FIM[1,3]

	# FIM[nu, p]
	FIM[2,3] = -1/(p*(z-1) + 1) * dlogzdnu
	FIM[3,2] = FIM[2,3]

	return(FIM)
}

# Compute the ZICMP information matrix using Monte Carlo.
fim.zicmp.mc = function(lambda, nu, p, reps)
{
	stopifnot(length(lambda) == 1)
	stopifnot(length(nu) == 1)
	stopifnot(length(p) == 1)
	grad.eps = getOption("COMPoissonReg.grad.eps", default = 1e-5)

	f = function(theta, data) {
		sum(dzicmp(data, lambda = theta[1], nu = theta[2], p = theta[3], log = TRUE))
	}

	FIM = matrix(0, 3, 3)
	x = rzicmp(reps, lambda, nu, p)
	for (r in 1:reps) {
		S = grad_fwd(f, x = c(lambda, nu, p), h = grad.eps, data = x[r])
		FIM = FIM + S %*% t(S)
	}

	return(FIM / reps)
}

# Compute the ZICMP information matrix when parameterized by regressions.
# If reps is NULL, attempt to use the exact expression (with numerical
# derivatives). Otherwise, use Monte Carlo approximation based on that
# many reps.
fim.zicmp.reg = function(X, S, W, beta, gamma, zeta, off.x, off.s, off.w, reps = NULL)
{
	n = nrow(X)
	qq = ncol(X) + ncol(S) + ncol(W)
	lambda = as.numeric(exp(X %*% beta + off.x))
	nu = as.numeric(exp(S %*% gamma + off.s))
	p = as.numeric(plogis(W %*% zeta + off.w))

	FIM = matrix(0, qq, qq)
	colnames(FIM) = c(
		sprintf("beta%d", 1:length(beta)),
		sprintf("gamma%d", 1:length(gamma)),
		sprintf("zeta%d", 1:length(zeta)))
	rownames(FIM) = colnames(FIM)

	# Compute the FIM with respect to each (lambda, nu, p), then
	# transform to the FIM of (beta, gamma, zeta)
	FIM.one.list = list()
	for (i in 1:n) {
		if (is.null(reps)) {
			FIM.one.list[[i]] = fim.zicmp(lambda[i], nu[i], p[i])
		} else {
			FIM.one.list[[i]] = fim.zicmp.mc(lambda[i], nu[i], p[i], reps)
		}
	}

	idx.beta = 1:length(beta)
	idx.gamma = 1:length(gamma) + length(beta)
	idx.zeta = 1:length(zeta) + length(gamma) + length(beta)

	# FIM[beta, beta]
	D = unlist(Map(function(x) { x[1,1] }, FIM.one.list))
	FIM[idx.beta, idx.beta] = t(X) %*% ((D * lambda^2) * X)

	# FIM[gamma, gamma]
	D = unlist(Map(function(x) { x[2,2] }, FIM.one.list))
	FIM[idx.gamma, idx.gamma] = t(S) %*% ((D * nu^2) * S)

	# FIM[zeta, zeta]
	D = unlist(Map(function(x) { x[3,3] }, FIM.one.list))
	FIM[idx.zeta, idx.zeta] = t(W) %*% ((D * p^2*(1-p)^2) * W)

	# FIM[beta, gamma]
	D = unlist(Map(function(x) { x[1,2] }, FIM.one.list))
	FIM[idx.beta, idx.gamma] = t(X) %*% (D * lambda * nu * S)
	FIM[idx.gamma, idx.beta] = t(FIM[idx.beta, idx.gamma])

	# FIM[zeta, gamma]
	D = unlist(Map(function(x) { x[3,2] }, FIM.one.list))
	FIM[idx.zeta, idx.gamma] = t(W) %*% (D * p*(1-p) * nu * S)
	FIM[idx.gamma, idx.zeta] = t(FIM[idx.zeta, idx.gamma])

	# FIM[beta, zeta]
	D = unlist(Map(function(x) { x[1,3] }, FIM.one.list))
	FIM[idx.zeta, idx.beta] = t(W) %*% (D * lambda * p*(1-p) * X)
	FIM[idx.beta, idx.zeta] = t(FIM[idx.zeta, idx.beta])

	return(FIM)
}
