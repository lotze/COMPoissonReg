fit.zicmp.reg = function(y, X, S, W, init, offset, fixed, control)
{
	start = Sys.time()
	d1 = ncol(X)
	d2 = ncol(S)
	d3 = ncol(W)
	n = length(y)
	qq = d1 + d2 + d3

	stopifnot("COMPoissonReg.init" %in% class(init))
	stopifnot("COMPoissonReg.offset" %in% class(offset))
	stopifnot("COMPoissonReg.fixed" %in% class(fixed))
	stopifnot("COMPoissonReg.control" %in% class(control))

	# Make sure dimensions match up
	stopifnot(n == nrow(X))
	stopifnot(n == nrow(S))
	stopifnot(n == nrow(W))
	stopifnot(n == length(offset$x))
	stopifnot(n == length(offset$s))
	stopifnot(n == length(offset$w))

	# Make sure fixed indices are between 1 and the corresponding dimension
	stopifnot(all(fixed$beta %in% seq_len(d1)))
	stopifnot(all(fixed$gamma %in% seq_len(d2)))
	stopifnot(all(fixed$zeta %in% seq_len(d3)))

	# Make sure initial values have the correct dimension
	if (d1 != length(init$beta)) {
		msg = sprintf("Length %d of init$beta is not equal to length %d of beta",
			length(init$beta), d1)
		stop(msg)
	}
	if (d2 != length(init$gamma)) {
		msg = sprintf("Length %d of init$gamma is not equal to length %d of gamma",
			length(init$gamma), d2)
		stop(msg)
	}
	if (d3 != length(init$zeta)) {
		msg = sprintf("Length %d of init$zeta is not equal to length %d of zeta",
			length(init$zeta), d3)
		stop(msg)
	}

	# Get settings from control
	optim.method = control$optim.method
	optim.control = control$optim.control
	hybrid.tol = control$hybrid.tol
	truncate.tol = control$truncate.tol
	ymax = control$ymax

	# "par" represents the vector of optimization variables; it contains only
	# coefficients which have not been fixed.
	#
	# "theta" is a list whose three elements represent the three vectors of
	# coefficients; these include fixed variables in addition to optimization
	# variables. Fixed variables are assumed to be set to their initial values.

	unfixed = get.fixed(
		beta = setdiff(seq_len(d1), fixed$beta),
		gamma = setdiff(seq_len(d2), fixed$gamma),
		zeta = setdiff(seq_len(d3), fixed$zeta)
	)

	idx.par1 = seq_along(unfixed$beta)
	idx.par2 = seq_along(unfixed$gamma) + length(unfixed$beta)
	idx.par3 = seq_along(unfixed$zeta) + length(unfixed$beta) +  length(unfixed$gamma)

	# The following functions transform back and forth between the "par" and
	# "theta" representations.
	
	par2theta = function(par) {
		beta = rep(NA, d1)
		beta[fixed$beta] = init$beta[fixed$beta]
		beta[unfixed$beta] = par[idx.par1]

		gamma = rep(NA, d2)
		gamma[fixed$gamma] = init$gamma[fixed$gamma]
		gamma[unfixed$gamma] = par[idx.par2]

		zeta = rep(NA, d3)
		zeta[fixed$zeta] = init$zeta[fixed$zeta]
		zeta[unfixed$zeta] = par[idx.par3]

		list(beta = beta, gamma = gamma, zeta = zeta)
	}

	theta2par = function(theta) {
		c(theta$beta[unfixed$beta],
			theta$gamma[unfixed$gamma],
			theta$zeta[unfixed$zeta])
	}

	loglik = function(par) {
		theta = par2theta(par)
		out = fitted.zicmp.internal(X, S, W, theta$beta, theta$gamma,
			theta$zeta, offset$x, offset$s, offset$w)
		loglik_zicmp(y, out$lambda, out$nu, out$p, hybrid.tol, truncate.tol, ymax)
	}

	if (!is.null(optim.control$fnscale)) {
		warning("optim.control$fnscale disregarded and taken as -1")
	}
	optim.control$fnscale = -1
	par.init = theta2par(init)
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
		sprintf("X:%s", colnames(X)[unfixed$beta]),
		sprintf("S:%s", colnames(S)[unfixed$gamma]),
		sprintf("W:%s", colnames(W)[unfixed$zeta]))

	loglik = res$value
	elapsed.sec = as.numeric(Sys.time() - start, units = "secs")

	res = list(theta.hat = theta.hat, H = H, opt.res = res,
		control = control, elapsed.sec = elapsed.sec, loglik = loglik, n = n,
		fixed = fixed, unfixed = unfixed)
	return(res)
}

# This is just provided as a convenience to call fit.zicmp.reg with some dummy
# values.
fit.cmp.reg = function(y, X, S, init, offset, fixed, control)
{
	n = length(y)
	W = matrix(NA, n, 0)
	init$zeta = numeric(0)
	offset$w = numeric(n)
	fixed$zeta = integer(0)

	fit.zicmp.reg(y, X, S, W, init, offset, fixed, control)
}
