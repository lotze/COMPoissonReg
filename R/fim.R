# These functions are currently not exported to users. Numerical computation
# of the FIM seems to be potentially unstable, whether computing by
# derivatives of the normalizing function, or explicitly programming the
# infinite sum expressions and truncating them. Instead, we currently
# provide variances via the Hessian from the optimizer.

# Compute the CMP information matrix using the exact expression, except
# use numerical derivatives of z functions. This avoids numerical issues
# with infinite sums.
fim_cmp = function(lambda, nu)
{
	f = function(theta) {
		z_hybrid(theta[1], theta[2], take_log = TRUE)
	}
	hess_eps = getOption("COMPoissonReg.hess.eps", default = 1e-2)
	H = hess_fwd(f, c(lambda, nu), h = hess_eps)
	H[1,1] = H[1,1] + cmp_expected_value(lambda, nu) / lambda^2
	return(H)
}

# Compute the CMP information matrix using Monte Carlo. This avoids
# issues with numerical second derivatives.
fim_cmp_mc = function(lambda, nu, reps)
{
	stopifnot(length(lambda) == 1)
	stopifnot(length(nu) == 1)
	grad_eps = getOption("COMPoissonReg.grad.eps", default = 1e-5)

	f = function(theta, data) {
		sum(dcmp(data, lambda = theta[1], nu = theta[2], log = TRUE))
	}

	FIM = matrix(0, 2, 2)
	x = rcmp(reps, lambda, nu)
	for (r in 1:reps) {
		S = grad_fwd(f, x = c(lambda, nu), h = grad_eps, data = x[r])
		FIM = FIM + S %*% t(S)
	}

	return(FIM / reps)
}

# Compute the ZICMP information matrix using the exact expression, except
# use numerical derivatives of z functions.
fim_zicmp = function(lambda, nu, p)
{
	stopifnot(length(lambda) == 1)
	stopifnot(length(nu) == 1)
	stopifnot(length(p) == 1)

	z = z_hybrid(lambda, nu)
	meany = zicmp_expected_value(lambda, nu, p)
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
	grad_eps = getOption("COMPoissonReg.grad.eps", default = 1e-5)
	hess_eps = getOption("COMPoissonReg.hess.eps", default = 1e-2)
	Dlogz = grad_fwd(f, c(lambda, nu), h = grad_eps)
	Hlogz = hess_fwd(f, c(lambda, nu), h = hess_eps)

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
fim_zicmp_mc = function(lambda, nu, p, reps)
{
	stopifnot(length(lambda) == 1)
	stopifnot(length(nu) == 1)
	stopifnot(length(p) == 1)
	grad_eps = getOption("COMPoissonReg.grad.eps", default = 1e-5)

	f = function(theta, data) {
		sum(dzicmp(data, lambda = theta[1], nu = theta[2], p = theta[3], log = TRUE))
	}

	FIM = matrix(0, 3, 3)
	x = rzicmp(reps, lambda, nu, p)
	for (r in 1:reps) {
		S = grad_fwd(f, x = c(lambda, nu, p), h = grad_eps, data = x[r])
		FIM = FIM + S %*% t(S)
	}

	return(FIM / reps)
}

# Compute the ZICMP information matrix when parameterized by regressions.
# If reps is NULL, attempt to use the exact expression (with numerical
# derivatives). Otherwise, use Monte Carlo approximation based on that
# many reps.
fim_zicmp_reg = function(X, S, W, beta, gamma, zeta, off_x, off_s, off_w, reps = NULL)
{
	n = nrow(X)
	qq = ncol(X) + ncol(S) + ncol(W)
	lambda = as.numeric(exp(X %*% beta + off_x))
	nu = as.numeric(exp(S %*% gamma + off_s))
	p = as.numeric(plogis(W %*% zeta + off_w))

	FIM = matrix(0, qq, qq)
	colnames(FIM) = c(
		sprintf("beta%d", 1:length(beta)),
		sprintf("gamma%d", 1:length(gamma)),
		sprintf("zeta%d", 1:length(zeta)))
	rownames(FIM) = colnames(FIM)

	# Compute the FIM with respect to each (lambda, nu, p), then
	# transform to the FIM of (beta, gamma, zeta)
	FIM_one_list = list()
	for (i in 1:n) {
		if (is.null(reps)) {
			FIM_one_list[[i]] = fim_zicmp(lambda[i], nu[i], p[i])
		} else {
			FIM_one_list[[i]] = fim_zicmp_mc(lambda[i], nu[i], p[i], reps)
		}
	}

	idx_beta = 1:length(beta)
	idx_gamma = 1:length(gamma) + length(beta)
	idx_zeta = 1:length(zeta) + length(gamma) + length(beta)

	# FIM[beta, beta]
	D = unlist(Map(function(x) { x[1,1] }, FIM_one_list))
	FIM[idx_beta, idx_beta] = t(X) %*% ((D * lambda^2) * X)

	# FIM[gamma, gamma]
	D = unlist(Map(function(x) { x[2,2] }, FIM_one_list))
	FIM[idx_gamma, idx_gamma] = t(S) %*% ((D * nu^2) * S)

	# FIM[zeta, zeta]
	D = unlist(Map(function(x) { x[3,3] }, FIM_one_list))
	FIM[idx_zeta, idx_zeta] = t(W) %*% ((D * p^2*(1-p)^2) * W)

	# FIM[beta, gamma]
	D = unlist(Map(function(x) { x[1,2] }, FIM_one_list))
	FIM[idx_beta, idx_gamma] = t(X) %*% (D * lambda * nu * S)
	FIM[idx_gamma, idx_beta] = t(FIM[idx_beta, idx_gamma])

	# FIM[zeta, gamma]
	D = unlist(Map(function(x) { x[3,2] }, FIM_one_list))
	FIM[idx_zeta, idx_gamma] = t(W) %*% (D * p*(1-p) * nu * S)
	FIM[idx_gamma, idx_zeta] = t(FIM[idx_zeta, idx_gamma])

	# FIM[beta, zeta]
	D = unlist(Map(function(x) { x[1,3] }, FIM_one_list))
	FIM[idx_zeta, idx_beta] = t(W) %*% (D * lambda * p*(1-p) * X)
	FIM[idx_beta, idx_zeta] = t(FIM[idx_zeta, idx_beta])

	return(FIM)
}
