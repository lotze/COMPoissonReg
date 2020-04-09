#' Supporting Functions for COM-Poisson Regression
#' 
#' @param object object of type \code{cmp}.
#' @param x object of type \code{cmp}.
#' @param k Penalty per parameter to be used in AIC calculation.
#' @param newdata New covariates to be used for prediction.
#' @param type Type of residual to be computed.
#' @param reps Number of bootstrap repetitions.
#' @param report.period Report progress every \code{report.period} iterations.
#' @param ... other arguments, such as \code{subset} and \code{na.action}.
#' 
#' @name glm.cmp, CMP support
NULL

#' @name glm.cmp, CMP support
#' @export
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
		z.value = round(z.val, 4),
		p.value = sprintf("%0.4g", p.val)
	)
	rownames(DF) = c(sprintf("X:%s", colnames(object$X)),
		sprintf("S:%s", colnames(object$S)))

	# If X, S, or W are intercept only, compute results for non-regression parameters
	DF.lambda = NULL
	DF.nu = NULL

	if (is.intercept.only(object$X) && is.zero.matrix(object$off.x)) {
		est = exp(object$beta)
		J = c(exp(object$beta), rep(0, d2))
		se = sqrt(t(J) %*% V %*% J)

		DF.lambda = data.frame(
			Estimate = round(est, 4),
			SE = round(se, 4)
		)
		rownames(DF.lambda) = "lambda"
	}

	if (is.intercept.only(object$S) && is.zero.matrix(object$off.s)) {
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

fitted.cmp.internal = function(X, S, beta, gamma, off.x, off.s)
{
	list(
		lambda = as.numeric(exp(X %*% beta + off.x)),
		nu = as.numeric(exp(S %*% gamma + off.s))
	)
}

#' @name glm.cmp, CMP support
#' @export
print.cmp = function(x, ...)
{
	printf("CMP coefficients\n")
	s = summary.cmp(x)
	tt = equitest(x)
	print(s$DF)

	if (!is.null(s$DF.lambda) || !is.null(s$DF.nu)) {
		printf("--\n")
		printf("Transformed intercept-only parameters\n")
		print(rbind(s$DF.lambda, s$DF.nu))
	}
	printf("--\n")
	printf("Chi-squared test for equidispersion\n")
	printf("X^2 = %0.4f, df = 1, ", tt$teststat)
	printf("p-value = %0.4e\n", tt$pvalue)
	printf("--\n")
	printf("Elapsed: %s   ", format_difftime(s$elapsed.sec))
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

#' @name glm.cmp, CMP support
#' @export
logLik.cmp = function(object, ...)
{
	object$loglik
}

#' @name glm.cmp, CMP support
#' @export
AIC.cmp = function(object, ..., k=2)
{
	-2*object$loglik + 2*length(coef(object))
}

#' @name glm.cmp, CMP support
#' @export
BIC.cmp = function(object, ...)
{
	n = length(object$y)
	-2*object$loglik + log(n)*length(coef(object))
}

#' @name glm.cmp, CMP support
#' @export
coef.cmp = function(object, ...)
{
	c(object$beta, object$gamma)
}

#' @name glm.cmp, CMP support
#' @export
nu.cmp = function(object, ...)
{
	out = fitted.cmp.internal(object$X, object$S, object$beta, object$gamma,
		object$off.x, object$off.s)
	out$nu
}

#' @name glm.cmp, CMP support
#' @export
sdev.cmp = function(object, ...)
{
	sqrt(diag(vcov(object)))
}

#' @name glm.cmp, CMP support
#' @export
vcov.cmp = function(object, ...)
{
	# Compute the covariance via Hessian from optimizer
	-solve(object$H)
}

#' @name equitest
#' @export
equitest.cmp = function(object, ...)
{
	if ("equitest" %in% names(object)) {
		return(object$equitest)
	}

	y = object$y
	d1 = ncol(object$X)
	d2 = ncol(object$S)

	out = fitted.cmp.internal(object$X, object$S, object$beta, object$gamma,
		object$off.x, object$off.s)
	lambda.hat = out$lambda
	nu.hat = out$nu

	out0 = fitted.cmp.internal(object$X, object$S, object$beta.glm,
		gamma = numeric(d2), object$off.x, object$off.s)
	lambda0.hat = out0$lambda

	logz = z_hybrid(lambda.hat, nu.hat, take_log = TRUE)
	teststat = -2 * sum(y*log(lambda0.hat) - lgamma(y+1) - lambda0.hat -
		y*log(lambda.hat) + nu.hat*lgamma(y+1) + logz)
	pvalue = pchisq(teststat, df = 1, lower.tail = FALSE)
	list(teststat = teststat, pvalue = pvalue)
}

#' @name glm.cmp, CMP support
#' @export
leverage.cmp = function(object, ...)
{
	y = object$y

	# 1) to code the WW matrix  (diagonal matrix with Var(Y_i) )
	ww = weights(object$X, object$S, object$beta, object$gamma, object$off.x, object$off.s)
	WW = diag(ww)

	out = fitted.cmp.internal(object$X, object$S, object$beta, object$gamma,
		object$off.x, object$off.s)
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

#' @name glm.cmp, CMP support
#' @export
deviance.cmp = function(object, ...)
{
	# Compute the COM-Poisson deviances exactly
	y = object$y
	n = length(y)

	#### Compute optimal log likelihood value for given nu-hat value
	beta.init = object$beta
	d1 = length(beta.init)
	ll.y = numeric(n)

	for (i in 1:n) {
		# loglik for single observation
		logf = function(beta) {
			out = fitted.cmp.internal(object$X[i,], object$S[i,], beta, object$gamma,
				object$off.x[i], object$off.s[i])
			y[i]*log(out$lambda) - out$nu*lgamma(y[i] + 1) -
				z_hybrid(out$lambda, out$nu, take_log = TRUE)
		}

		# Determine the MLEs
		# Using optim rather than nlm because optim handles -Inf more gracefully, if encountered
		opt.method = getOption("COMPoissonReg.optim.method")
		res = optim(beta.init, logf, method = opt.method, control = list(fnscale = -1))
		ll.y[i] = res$value
	}

	# Compute exact deviances
	out = fitted.cmp.internal(object$X, object$S, object$beta, object$gamma,
		object$off.x, object$off.s)
	ll.mu = y*log(out$lambda) - out$nu * lgamma(y+1) -
		z_hybrid(out$lambda, out$nu, take_log = TRUE)
	d = -2*(ll.mu - ll.y)
	lev = leverage.cmp(object)
	cmpdev = d / sqrt(1 - lev)
	return(cmpdev)
}

#' @name glm.cmp, CMP support
#' @export
residuals.cmp = function(object, type = c("raw", "quantile"), ...)
{
	out = fitted.cmp.internal(object$X, object$S, object$beta, object$gamma, object$off.x, object$off.s)
	y.hat = ecmp(out$lambda, out$nu)

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

#' @name glm.cmp, CMP support
#' @export
predict.cmp = function(object, newdata = NULL, ...)
{
	if (is.null(newdata)) {
		X = object$X
		S = object$S
		off.x = object$off.x
		off.s = object$off.s
	} else {
		# Only attempt to process newdata as a data.frame
		newdata = as.data.frame(newdata)

		# If the response was not included in newdata, add a column with zeros
		response_name = all.vars(object$formula.lambda)[1]
		if (is.null(newdata[[response_name]])) {
			newdata[[response_name]] = 0
		}

		mf.x = model.frame(object$formula.lambda, data = newdata, ...)
		mf.s = model.frame(object$formula.nu, data = newdata, ...)
		X = model.matrix(object$formula.lambda, mf.x)
		S = model.matrix(object$formula.nu, mf.s)
		off.x = model.offset(mf.x)
		off.s = model.offset(mf.s)

		n.new = nrow(X)
		if (is.null(off.x)) { off.x = rep(0, n.new) }
		if (is.null(off.s)) { off.s = rep(0, n.new) }
	}

	out = fitted.cmp.internal(X, S, object$beta, object$gamma, off.x, off.s)
	ecmp(out$lambda, out$nu)
}

#' @name glm.cmp, CMP support
#' @export
parametric.bootstrap.cmp = function(object, reps = 1000, report.period = reps+1, ...)
{
	n = nrow(object$X)
	d1 = ncol(object$X)
	d2 = ncol(object$S)

	out = fitted.cmp.internal(object$X, object$S, object$beta, object$gamma, object$off.x, object$off.s)
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
				off.x = object$off.x, off.s = object$off.s)
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

weights = function(X, S, beta, gamma, off.x, off.s)
{
	out = fitted.cmp.internal(X, S, beta, gamma, off.x, off.s)

	# Compute the parts that comprise the weight functions
	w1 = z_prodj2(out$lambda, out$nu)
	w2 = z_prodj(out$lambda, out$nu)
	w3 = z_hybrid(out$lambda, out$nu)

	Ey2 = w1 / w3
	E2y = (w2 / w3)^2

	# Determine the weights
	weight = Ey2 - E2y
	return(weight)
}

constantCMPfitsandresids = function(lambda.hat, nu.hat, y = 0)
{
	# Determine estimated lambda.hat, fit, and residuals
	lambda.hat = as.numeric(lambda.hat)
	nu.hat = as.numeric(nu.hat)
	fit = lambda.hat^(1/nu.hat) - ((nu.hat - 1)/(2*nu.hat))
	resid = y - fit
	list(fit = fit, resid = resid)
}

cmp.mse = function(y, X, S, beta, gamma, off.x, off.s)
{
	out = fitted.cmp.internal(X, S, beta, gamma, off.x, off.s)
	fnr = constantCMPfitsandresids(out$lambda, out$nu, y)
	res = fnr$resid
	mse = mean(res^2)
	return(mse)
}
