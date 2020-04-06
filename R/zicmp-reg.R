#' Supporting Functions for ZICMP Regression
#' 
#' @param object object of type \code{zicmp}.
#' @param x object of type \code{zicmp}.
#' @param k Penalty per parameter to be used in AIC calculation.
#' @param newdata New covariates to be used for prediction.
#' @param type Type of residual to be computed.
#' @param reps Number of bootstrap repetitions.
#' @param report.period Report progress every \code{report.period} iterations.
#' @param ... other model parameters, such as data.
#' 
#' @name glm.cmp, ZICMP support
NULL

#' @name glm.cmp, ZICMP support
#' @export
summary.zicmp = function(object, ...)
{
	n = nrow(object$X)
	d1 = ncol(object$X)
	d2 = ncol(object$S)
	d3 = ncol(object$W)

	V = vcov(object)
	est = coef(object)
	se = sdev(object)
	z.val = est / se
	p.val = 2*pnorm(-abs(z.val))
	qq = length(est)

	DF = data.frame(
		Estimate = round(est, 4),
		SE = round(se, 4),
		z.value = round(z.val, 4),
		p.value = sprintf("%0.4g", p.val)
	)
	rownames(DF) = c(
		sprintf("X:%s", colnames(object$X)),
		sprintf("S:%s", colnames(object$S)),
		sprintf("W:%s", colnames(object$W))
	)

	# If X, S, or W are intercept only, compute results for non-regression parameters
	# Exclude offsets from these calculations
	DF.lambda = NULL
	DF.nu = NULL
	DF.p = NULL

	if (is.intercept.only(object$X) && is.zero.matrix(object$off.x)) {
		est = exp(object$beta)
		J = c(exp(object$beta), rep(0, d2), rep(0, d3))
		se = sqrt(t(J) %*% V %*% J)

		DF.lambda = data.frame(
			Estimate = round(est, 4),
			SE = round(se, 4)
		)
		rownames(DF.lambda) = "lambda"
	}

	if (is.intercept.only(object$S) && is.zero.matrix(object$off.s)) {
		est = exp(object$gamma)
		J = c(rep(0, d1), exp(object$gamma), rep(0, d3))
		se = sqrt(t(J) %*% V %*% J)

		DF.nu = data.frame(
			Estimate = round(est, 4),
			SE = round(se, 4)		)
		rownames(DF.nu) = "nu"
	}

	if (is.intercept.only(object$W) && is.zero.matrix(object$off.w)) {
		est = plogis(object$zeta)
		J = c(rep(0, d1), rep(0, d2), dlogis(object$zeta))
		se = sqrt(t(J) %*% V %*% J)

		DF.p = data.frame(
			Estimate = round(est, 4),
			SE = round(se, 4)
		)
		rownames(DF.p) = "p"
	}

	list(DF = DF, DF.lambda = DF.lambda, DF.nu = DF.nu, DF.p = DF.p,
		n = length(object$y),
		loglik = logLik(object),
		aic = AIC(object),
		bic = BIC(object),
		opt.method = object$opt.method,
		opt.res = object$opt.res,
		elapsed.sec = object$elapsed.sec
	)
}

fitted.zicmp.internal = function(X, S, W, beta, gamma, zeta, off.x, off.s, off.w)
{
	list(
		lambda = as.numeric(exp(X %*% beta + off.x)),
		nu = as.numeric(exp(S %*% gamma + off.s)),
		p = as.numeric(plogis(W %*% zeta + off.w))
	)
}

#' @name glm.cmp, ZICMP support
#' @export
print.zicmp = function(x, ...)
{
	printf("ZICMP coefficients\n")
	s = summary(x)
	tt = equitest(x)
	print(s$DF)

	if (!is.null(s$DF.lambda) || !is.null(s$DF.nu) || !is.null(s$DF.p)) {
		printf("--\n")
		printf("Transformed intercept-only parameters\n")
		print(rbind(s$DF.lambda, s$DF.nu, s$DF.p))
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
	printf("Converged status: %d   ", s$opt.convergence)
	printf("Message: %s\n", s$opt.message)
}

#' @name glm.cmp, ZICMP support
#' @export
logLik.zicmp = function(object, ...)
{
	object$loglik
}

#' @name glm.cmp, ZICMP support
#' @export
AIC.zicmp = function(object, ..., k = 2)
{
	-2*object$loglik + 2*length(coef(object))
}

#' @name glm.cmp, ZICMP support
#' @export
BIC.zicmp = function(object, ...)
{
	n = length(object$y)
	-2*object$loglik + log(n)*length(coef(object))
}

#' @name glm.cmp, ZICMP support
#' @export
coef.zicmp = function(object, ...)
{
	c(object$beta, object$gamma, object$zeta)
}

#' @name glm.cmp, ZICMP support
#' @export
nu.zicmp = function(object, ...)
{
	out = fitted.zicmp.internal(object$X, object$S, object$W, object$beta,
		object$gamma, object$zeta, object$off.x, object$off.s, object$off.w)
	out$nu
}

#' @name glm.cmp, ZICMP support
#' @export
sdev.zicmp = function(object, ...)
{
	sqrt(diag(vcov(object)))
}

#' @name glm.cmp, ZICMP support
#' @export
vcov.zicmp = function(object, ...)
{
	# Compute the covariance via Hessian from optimizer
	-solve(object$H)
}

#' @name equitest
#' @export
equitest.zicmp = function(object, ...)
{
	if ("equitest" %in% names(object)) {
		return(object$equitest)
	}

	y = object$y
	X = object$X
	S = object$S
	W = object$W
	beta.hat = object$beta
	gamma.hat = object$gamma
	zeta.hat = object$zeta
	off.x = object$off.x
	off.s = object$off.s
	off.w = object$off.w
	n = length(y)
	d1 = ncol(X)
	d2 = ncol(S)
	d3 = ncol(W)

	fit0.out = fit.zip.reg(y, X, W, beta.init = beta.hat,
		zeta.init = zeta.hat, off.x = off.x, off.w = off.w)
	beta0.hat = fit0.out$theta.hat$beta
	zeta0.hat = fit0.out$theta.hat$zeta

	out1 = fitted.zicmp.internal(X, S, W, beta.hat, gamma.hat, zeta.hat, off.x, off.s, off.w)
	lambda.hat = out1$lambda
	nu.hat = out1$nu
	p.hat = out1$p

	out0 = fitted.zicmp.internal(X, S, W, beta0.hat, numeric(d2), zeta0.hat, off.x, off.s, off.w)
	lambda0.hat = out0$lambda
	p0.hat = out0$p

	logff = dzicmp(y, lambda.hat, nu.hat, p.hat, log = TRUE)
	logff0 = dzip(y, lambda0.hat, p0.hat, log = TRUE)

	X2 = 2*(sum(logff) - sum(logff0))
	pvalue = pchisq(X2, df = d3, lower.tail = FALSE)
	list(teststat = X2, pvalue = pvalue)
}

#' @name glm.cmp, ZICMP support
#' @export
leverage.zicmp = function(object, ...)
{
	stop("This function is not yet implemented")
}

#' @name glm.cmp, ZICMP support
#' @export
deviance.zicmp = function(object, ...)
{
	y = object$y
	X = object$X
	S = object$S
	W = object$W
	n = length(y)
	d1 = ncol(X)
	d2 = ncol(S)
	d3 = ncol(W)
	off.x = object$off.x
	off.s = object$off.s
	off.w = object$off.w

	par.hat = c(object$beta, object$gamma, object$zeta)
	par.init = par.hat
	ll = numeric(n)
	ll.star = numeric(n)

	for (i in 1:n) {
		loglik = function(par){
			beta = par[1:d1]
			gamma = par[1:d2 + d1]
			zeta = par[1:d3 + d1 + d2]
			out = fitted.zicmp.internal(X[i,], S[i,], W[i,], beta, gamma, zeta,
				off.x[i], off.s[i], off.w[i])
			dzicmp(y[i], out$lambda, out$nu, out$p, log = TRUE)
		}

		# Maximize loglik for ith obs
		res = optim(par.init, loglik, control = list(fnscale = -1))
		ll.star[i] = res$value

		# loglik maximized over all the obs, evaluated for ith obs
		ll[i] = loglik(par.hat)
	}

	dev = -2*(ll - ll.star)
	leverage = leverage(object)
	cmpdev = dev / sqrt(1 - leverage)
	return(cmpdev)
}

#' @name glm.cmp, ZICMP support
#' @export
residuals.zicmp = function(object, type = c("raw", "quantile"), ...)
{
	out = fitted.zicmp.internal(object$X, object$S, object$W, object$beta,
		object$gamma, object$zeta, object$off.x, object$off.s, object$off.w)
	lambda.hat = out$lambda
	nu.hat = out$nu
	p.hat = out$p
	y.hat = predict.zicmp(object)

	type = match.arg(type)
	if (type == "raw") {
		res = object$y - y.hat
	} else if (type == "quantile") {
		res = rqres.zicmp(object$y, lambda.hat, nu.hat, p.hat)
	} else {
		stop("Unsupported residual type")
	}

	return(as.numeric(res))
}

#' @name glm.cmp, ZICMP support
#' @export
predict.zicmp = function(object, newdata = NULL, ...)
{
	if (is.null(newdata)) {
		X = object$X
		S = object$S
		W = object$W
		off.x = object$off.x
		off.s = object$off.s
		off.w = object$off.w
	} else {
		# If any of the original models had an intercept added via model.matrix, they
		# will have an "(Intercept)" column. Let's add an "(Intercept)" to newdata
		# in case the user didn't make one.
		newdata$'(Intercept)' = 1

		n.new = nrow(newdata)
		X = as.matrix(newdata[,colnames(object$X)])
		S = as.matrix(newdata[,colnames(object$S)])
		W = as.matrix(newdata[,colnames(object$W)])

		mf.x = model.frame(object$formula.lambda, data = newdata, ...)
		mf.s = model.frame(object$formula.nu, data = newdata, ...)
		mf.p = model.frame(object$formula.p, data = newdata, ...)
		off.x = model.offset(mf.x)
		off.s = model.offset(mf.s)
		off.w = model.offset(mf.p)
		if (is.null(off.x)) { off.x = rep(0, n.new) }
		if (is.null(off.s)) { off.s = rep(0, n.new) }
		if (is.null(off.w)) { off.w = rep(0, n.new) }
	}

	out = fitted.zicmp.internal(X, S, W, object$beta,
		object$gamma, object$zeta, off.x, off.s, off.w)
	y.hat = zicmp.expected.value(out$lambda, out$nu, out$p)
	return(y.hat)
}

#' @name glm.cmp, ZICMP support
#' @export
parametric.bootstrap.zicmp = function(object, reps = 1000, report.period = reps+1, ...)
{
	n = length(object$y)
	qq = length(object$beta) + length(object$gamma) + length(object$zeta)
	theta.boot = matrix(NA, reps, qq)

	out = fitted.zicmp.internal(object$X, object$S, object$W, object$beta,
		object$gamma, object$zeta, object$off.x, object$off.s, object$off.w)
	lambda.hat = out$lambda
	nu.hat = out$nu
	p.hat = out$p

	for (r in 1:reps) {
		if (r %% report.period == 0) {
			logger("Starting bootstrap rep %d\n", r)
		}

		# Generate bootstrap samples of the full dataset using MLE
		y.boot = rzicmp(n, lambda.hat, nu.hat, p.hat)

		# Take each of the bootstrap samples and fit model to generate bootstrap
		# estimates
		tryCatch({
			fit.boot = fit.zicmp.reg(y.boot, object$X, object$S, object$W,
				object$beta.init, object$gamma.init, object$zeta.init,
				off.x = object$off.x, off.s = object$off.s, off.w = object$off.w)
			theta.boot[r,] = unlist(fit.boot$theta.hat)
		},
		error = function(e) {
			# Do nothing now; emit a warning later
		})
	}

	cnt = sum(rowSums(is.na(theta.boot)) > 0)
	if (cnt > 0) {
		warning(sprintf("%d out of %d bootstrap iterations failed", cnt, reps))
	}

	colnames(theta.boot) = c(
		sprintf("X:%s", colnames(object$X)),
		sprintf("S:%s", colnames(object$S)),
		sprintf("W:%s", colnames(object$W))
	)

	return(theta.boot)
}
