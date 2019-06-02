summary.cmp = function(object, ...)
{
	n = nrow(object$X)
	d1 = ncol(object$X)
	d2 = ncol(object$S)

	V = vcov(object)
	est = coef(object)
	se = sdev(object)
	z.val = est / se
	p.val = 2*pnorm(-abs(z.val))

	DF = data.frame(
		Estimate = round(est, 4),
		SE = round(se, 4),
		z_value = round(z.val, 4),
		p_value = sprintf("%0.4g", p.val)
	)
	rownames(DF) = c(sprintf("X:%s", colnames(object$X)),
		sprintf("S:%s", colnames(object$S)))

	# If X, S, or W are intercept only, compute results for non-regression parameters
	DF.lambda = NULL
	DF.nu = NULL
	X.int = matrix(1, n, 1)

	is.intercept.only = function(A, eps = 1e-12) {
		prod(dim(A) == c(n,1)) & norm(A - 1, type = "F") < eps
	}

	if (is.intercept.only(object$X)) {
		est = exp(object$beta)
		J = c(exp(object$beta), rep(0, d2))
		se = sqrt(t(J) %*% V %*% J)

		DF.lambda = data.frame(
			Estimate = round(est, 4),
			SE = round(se, 4)
		)
		rownames(DF.lambda) = "lambda"
	}

	if (is.intercept.only(object$S)) {
		est = exp(object$gamma)
		J = c(rep(0, d1), exp(object$gamma))
		se = sqrt(t(J) %*% V %*% J)

		DF.nu = data.frame(
			Estimate = round(est, 4),
			SE = round(se, 4)
		)
		rownames(DF.nu) = "nu"
	}

	list(DF = DF, DF.lambda = DF.lambda, DF.nu = DF.nu,
		n = n,
		loglik = logLik(object),
		aic = AIC(object),
		bic = BIC(object),
		opt.method = object$opt.method,
		opt.message = object$opt.res$message,
		opt.convergence = object$opt.res$convergence,
		elapsed.sec = object$elapsed.sec
	)
}

fitted_cmp_internal = function(X, S, beta, gamma, off.X, off.S)
{
	list(
		lambda = as.numeric(exp(X %*% beta + off.X)),
		nu = as.numeric(exp(S %*% gamma + off.S))
	)
}

print.cmp = function(x, ...)
{
	printf("CMP coefficients\n")
	s = summary.cmp(x)
	tt = equitest(x)
	print(s$DF)

	if (!is.null(s$DF.lambda) || !is.null(s$DF.nu)) {
		printf("--\n")
		printf("Transformed intercept-only parameters (excluding offsets)\n")
		print(rbind(s$DF.lambda, s$DF.nu))
	}
	printf("--\n")
	printf("Chi-squared test for equidispersion\n")
	printf("X^2 = %0.4f, df = 1, ", tt$teststat)
	printf("p-value = %0.4e\n", tt$pvalue)
	printf("--\n")
	printf("Elapsed Sec: %0.2f   ", s$elapsed.sec)
	printf("Sample size: %d   ", s$n)
	printf("SEs via Hessian\n")
	printf("LogLik: %0.4f   ", s$loglik)
	printf("AIC: %0.4f   ", s$aic)
	printf("BIC: %0.4f   ", s$bic)
	printf("\n")
	printf("Optimization Method: %s   ", s$opt.method)
	printf("Converged status: %d\n", s$opt.convergence)
	printf("Message: %s\n", s$opt.message)
}

logLik.cmp = function(object, ...)
{
	object$loglik
}

AIC.cmp = function(object, ..., k=2)
{
	-2*object$loglik + 2*length(coef(object))
}

BIC.cmp = function(object, ...)
{
	n = length(object$y)
	-2*object$loglik + log(n)*length(coef(object))
}

coef.cmp = function(object, ...)
{
	c(object$beta, object$gamma)
}

nu.cmp = function(object, ...)
{
	out = fitted_cmp_internal(object$X, object$S, object$beta, object$gamma,
		object$off.X, object$off.S)
	out$nu
}

sdev.cmp = function(object, ...)
{
	sqrt(diag(vcov(object)))
}

vcov.cmp = function(object, ...)
{
	# Compute the covariance via Hessian from optimizer
	-solve(object$H)
}

equitest.cmp = function(object, ...)
{
	if ("equitest" %in% names(object)) {
		return(object$equitest)
	}

	y = object$y
	d1 = ncol(object$X)
	d2 = ncol(object$S)

	out = fitted_cmp_internal(object$X, object$S, object$beta, object$gamma,
		object$off.X, object$off.S)
	lambda.hat = out$lambda
	nu.hat = out$nu
	
	out0 = fitted_cmp_internal(object$X, object$S, object$beta.glm,
		gamma = numeric(d2), object$off.X, object$off.S)
	lambda0.hat = out0$lambda

	logz = z_hybrid(lambda.hat, nu.hat, take_log = TRUE)
	teststat = -2 * sum(y*log(lambda0.hat) - lgamma(y+1) - lambda0.hat -
		y*log(lambda.hat) + nu.hat*lgamma(y+1) + logz)
	pvalue = pchisq(teststat, df = 1, lower.tail = FALSE)
	list(teststat = teststat, pvalue = pvalue)
}

leverage.cmp = function(object, ...)
{
	y = object$y

	# 1) to code the WW matrix  (diagonal matrix with Var(Y_i) )
	ww = weights(object$X, object$S, object$beta, object$gamma, object$off.X, object$off.S)
	WW = diag(ww)

	out = fitted_cmp_internal(object$X, object$S, object$beta, object$gamma,
		object$off.X, object$off.S)
	lambda.hat = out$lambda
	nu.hat = out$nu

	#    and X matrix (in Appendix)
	E.y = z_prodj(lambda.hat, nu.hat) / z_hybrid(lambda.hat, nu.hat)
	E.logfacty = z_prodlogj(lambda.hat, nu.hat) / z_hybrid(lambda.hat, nu.hat)
	extravec = (-lgamma(y+1) + E.logfacty)/(y - E.y)
	curlyX.mat = cbind(object$X, extravec)

	# 2) to compute H using eq (12)  on p. 11
	H1 = t(curlyX.mat) %*% sqrt(WW)
	H2 = solve(t(curlyX.mat) %*% WW %*% curlyX.mat)
	H = t(H1) %*% H2 %*% H1
	diagH = diag(H)
	return(diagH)
}

deviance.cmp = function(object, ...)
{
	# Compute the COM-Poisson deviances exactly
	y = object$y
	n = length(y)
	leverage = leverage.cmp(object)

	#### Compute optimal log likelihood value for given nu-hat value
	beta.init = object$beta
	d1 = length(beta.init)
	OptimalLogLi = rep(0,n)
	iterct = rep(0,n)

	for (i in 1:length(y)){
		# Create -logL = -logf (because considering single observation) so that we
		# take the minimum of this function (which equals the max of logL)
		logf = function(par){
			beta = par[1:d1]
			out = fitted_cmp_internal(object$X[i,], object$S[i,], beta, object$gamma,
				object$off.X[i], object$off.S[i])
			y[i]*log(out$lambda) - out$nu*lgamma(y[i]+1) -
				z_hybrid(out$lambda, out$nu, take_log = TRUE)
		}

		# Determine the MLEs
		# Using optim rather than nlm because optim handles -Inf more gracefully, if encountered
		BetaEstResult = optim(beta.init, logf, method = "L-BFGS-B", control = list(fnscale = -1))
		OptimalLogLi[i] = BetaEstResult$value
	}

	#### Compute exact deviances
	out = fitted_cmp_internal(object$X, object$S, object$beta, object$gamma, object$off.X, object$off.S)
	OptimalLogL.mu = (y*log(out$lambda)) - (out$nu * lgamma(y+1)) -
		z_hybrid(out$lambda, out$nu, take_log = TRUE)
	OptimalLogL.y = OptimalLogLi
	d = -2*(OptimalLogL.mu - OptimalLogL.y)
	cmpdev = d/(sqrt(1-leverage))
	return(cmpdev)
}

residuals.cmp = function(object, type = c("raw", "quantile"), ...)
{
	out = fitted_cmp_internal(object$X, object$S, object$beta, object$gamma, object$off.X, object$off.S)
	y.hat = predict.cmp(object, newdata = object$X)

	type = match.arg(type)
	if (type == "raw") {
		res = object$y - y.hat
	} else if (type == "quantile") {
		res = rqres.cmp(object$y, lambda = out$lambda, nu = out$nu)
	} else {
		stop("Unsupported residual type")
	}

	return(as.numeric(res))
}

predict.cmp = function(object, newdata = NULL, ...)
{
	if (!is.null(newdata)) {
		# If any of the original models had an intercept added via model.matrix, they
		# will have an "(Intercept)" column. Let's add an "(Intercept)" to newdata
		# in case the user didn't make one.
		newdata = as.data.frame(newdata)
        newdata$'(Intercept)' = 1
		X = as.matrix(newdata[,colnames(object$X)])
		S = as.matrix(newdata[,colnames(object$S)])
	} else {
		X = object$X
		S = object$S
	}

	out = fitted_cmp_internal(X, S, object$beta, object$gamma, object$off.X, object$off.S)
	fnr = constantCMPfitsandresids(out$lambda, out$nu)
	y.hat = fnr$fit
	return(y.hat)
}

parametric_bootstrap.cmp = function(object, reps = 1000, report.period = reps+1, ...)
{
	n = nrow(object$X)
	d1 = ncol(object$X)
	d2 = ncol(object$S)

	out = fitted_cmp_internal(object$X, object$S, object$beta, object$gamma, object$off.X, object$off.S)
	lambda.hat = out$lambda
	nu.hat = out$nu

	# Generate `reps` samples, using beta.hat and nu.hat from full dataset
	# Run CMP regression on each boostrap sample to generate new beta and nu estimates
	
	boot.out = matrix(NA, nrow = reps, ncol = d1 + d2)

	for (r in 1:reps){
		if (r %% report.period == 0) {
			logger("Starting bootstrap rep %d\n", r)
		}
		y.boot = rcmp(n, lambda.hat, nu.hat)
		tryCatch({
			res = fit.cmp.reg(y = y.boot, X = object$X, S = object$S,
				beta.init = object$beta.glm, gamma.init = object$gamma,
				off.X = object$off.X, off.S = object$off.S)
			boot.out[r,] = unlist(res$theta.hat)
		}, error = function(e) {
			# Do nothing now; emit a warning later
		})
	}

	cnt = sum(rowSums(is.na(boot.out)) > 0)
	if (cnt > 0) {
		warning(sprintf("%d out of %d bootstrap iterations failed", cnt, reps))
	}

	colnames(boot.out) = c(colnames(object$X), colnames(object$S), recursive=TRUE)
	return(boot.out)
}

constantCMPfitsandresids = function(lambda.hat, nu.hat, y=0)
{
	# Determine estimated lambda.hat, fit, and residuals
	lambda.hat = as.numeric(lambda.hat)
	nu.hat = as.numeric(nu.hat)
	fit = lambda.hat^(1/nu.hat) - ((nu.hat - 1)/(2*nu.hat))
	resid = y - fit
	list(fit = fit, resid = resid)
}

weights = function(X, S, beta, gamma, off.X, off.S)
{
	out = fitted_cmp_internal(X, S, beta, gamma, off.X, off.S)
	lambda = out$lambda
	nu = out$nu

	# Compute the parts that comprise the weight functions
	w1 = z_prodj2(lambda, nu)
	w2 = z_prodj(lambda, nu)
	w3 = z_hybrid(lambda, nu)

	Ey2 = w1/w3
	E2y = (w2/w3)^2

	# Determine the weights
	weight = Ey2 - E2y
	return(weight)
}

CMP.MSE = function(y, X, S, beta, gamma, off.X, off.S)
{
	out = fitted_cmp_internal(X, S, beta, gamma, off.X, off.S)
	lambda = out$lambda
	nu = out$nu

	fnr = constantCMPfitsandresids(lambda, nu, y)
	res = fnr$resid
	MSE = mean(res^2)
	return(MSE)
}
