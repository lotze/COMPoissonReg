#' Supporting Functions for ZICMP Regression
#' 
#' @param object object of type \code{zicmp}.
#' @param x object of type \code{zicmp}.
#' @param k Penalty per parameter to be used in AIC calculation.
#' @param newdata New covariates to be used for prediction.
#' @param type Specifies quantity to be computed. See details.
#' @param reps Number of bootstrap repetitions.
#' @param report.period Report progress every \code{report.period} iterations.
#' @param ... other arguments, such as \code{subset} and \code{na.action}.
#' 
#' @details
#' The function \code{residuals} returns raw residuals when
#' \code{type = "raw"}  and quantile residuals when
#' \code{type = "quantile"}.
#' 
#' The function \code{predict} returns expected values of the outcomes,
#' eveluated at the computed estimates, when \code{type = "response"}. When
#' \code{type = "link"}, a \code{data.frame} is instead returned with
#' columns corresponding to estimates of \code{lambda}, \code{nu}, and
#' \code{p}.
#' 
#' The function \code{coef} returns a vector of coefficient estimates in
#' the form \code{c(beta, gamma, zeta)} when \code{type = "vector"}. When
#' \code{type = "list"}, the estimates are returned as a list with named
#' elements \code{beta} and \code{gamma}, and \code{zeta}.
#' 
#' The \code{type} argument behaves the same for the \code{sdev} function
#' as it does for \code{coef}.
#' 
#' @name glm.cmp, ZICMP support
NULL

#' @name glm.cmp, ZICMP support
#' @export
summary.zicmpfit = function(object, ...)
{
	n = nrow(object$X)
	d1 = ncol(object$X)
	d2 = ncol(object$S)
	d3 = ncol(object$W)
	qq = d1 + d2 + d3

	# We need the indices of the fixed coefficients and the ones included in
	# optimization.
	fixed = object$fixed
	unfixed = object$unfixed
	idx.par1 = seq_along(unfixed$beta)
	idx.par2 = seq_along(unfixed$gamma) + length(unfixed$beta)
	idx.par3 = seq_along(unfixed$zeta) + length(unfixed$beta) +  length(unfixed$gamma)

	V = vcov(object)
	est = coef(object)

	# In the vector of SEs, include an NA entry if the variable was fixed
	# This NA will propagate to the corresponding z-value and p-value as well.
	se = rep(NA, qq)
	se[unfixed$beta] = sdev(object)[idx.par1]
	se[unfixed$gamma + d1] = sdev(object)[idx.par2]
	se[unfixed$zeta + d1 + d2] = sdev(object)[idx.par3]

	z.val = est / se
	p.val = 2*pnorm(-abs(z.val))
	qq = length(est)

	DF = data.frame(
		Estimate = round(est, 4),
		SE = round(se, 4),
		z.value = round(z.val, 4),
		p.value = sprintf("%0.4g", p.val)
	)
	colnames(DF) = c("Estimate", "SE", "z-value", "p-value")
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

	# In each block below, make sure to consider only the non-fixed variables
	# for the Jacobian and Hessian. If one of the intercepts was fixed, it should
	# result in an SE of zero.
	
	if (is.intercept.only(object$X) && is.zero.matrix(object$offset$x)) {
		if (length(fixed$beta) > 0) {
			se = 0
		} else {
			J = c(exp(object$beta), numeric(length(idx.par2)), numeric(length(idx.par3)))
			se = sqrt(t(J) %*% V %*% J)
		}
		est = exp(object$beta)		
		DF.lambda = data.frame(
			Estimate = round(est, 4),
			SE = round(se, 4)
		)
		rownames(DF.lambda) = "lambda"
	}

	if (is.intercept.only(object$S) && is.zero.matrix(object$offset$s)) {
		if (length(fixed$gamma) > 0) {
			se = 0
		} else {
			J = c(numeric(length(idx.par1)), exp(object$gamma), numeric(rep(length(idx.par3))))
			se = sqrt(t(J) %*% V %*% J)
		}
		est = exp(object$gamma)
		DF.nu = data.frame(
			Estimate = round(est, 4),
			SE = round(se, 4)		)
		rownames(DF.nu) = "nu"
	}

	if (is.intercept.only(object$W) && is.zero.matrix(object$offset$w)) {
		if (length(fixed$zeta) > 0) {
			se = 0
		} else {
			J = c(numeric(length(idx.par1)), numeric(length(idx.par2)), dlogis(object$zeta))
			se = sqrt(t(J) %*% V %*% J)
		}
		est = plogis(object$zeta)
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
		optim.method = object$control$optim.method,
		opt.message = object$opt.res$message,
		opt.convergence = object$opt.res$convergence,
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
print.zicmpfit = function(x, ...)
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
	if (is.character(tt)) {
		printf("--\n")
		cat(paste(tt, collapse = "\n"))
		printf("\n")
	} else {
		printf("--\n")
		printf("Chi-squared test for equidispersion\n")
		printf("X^2 = %0.4f, df = %d, ", tt$teststat, tt$df)
		printf("p-value = %0.4e\n", tt$pvalue)
	}
	printf("--\n")
	printf("Elapsed: %s   ", format.difftime(s$elapsed.sec))
	printf("Sample size: %d   ", s$n)
	printf("%s interface\n", x$interface)
	printf("LogLik: %0.4f   ", s$loglik)
	printf("AIC: %0.4f   ", s$aic)
	printf("BIC: %0.4f   ", s$bic)
	printf("\n")
	printf("Optimization Method: %s   ", s$optim.method)
	printf("Converged status: %d\n", s$opt.convergence)
	printf("Message: %s\n", s$opt.message)
}

#' @name glm.cmp, ZICMP support
#' @export
logLik.zicmpfit = function(object, ...)
{
	object$loglik
}

#' @name glm.cmp, ZICMP support
#' @export
AIC.zicmpfit = function(object, ..., k = 2)
{
	-2*object$loglik + 2*length(coef(object))
}

#' @name glm.cmp, ZICMP support
#' @export
BIC.zicmpfit = function(object, ...)
{
	n = length(object$y)
	-2*object$loglik + log(n)*length(coef(object))
}

#' @name glm.cmp, ZICMP support
#' @export
coef.zicmpfit = function(object, type = c("vector", "list"), ...)
{
	switch(match.arg(type),
		vector = c(object$beta, object$gamma, object$zeta),
		list = list(beta = object$beta, gamma = object$gamma,
			zeta = object$zeta)
	)
}

#' @name glm.cmp, ZICMP support
#' @export
nu.zicmpfit = function(object, ...)
{
	# This function is deprecated - use predict instead
	.Deprecated("predict(object, type = \"link\")")
	link = predict(object, type = "link")
	link$nu
}

#' @name glm.cmp, ZICMP support
#' @export
sdev.zicmpfit = function(object, type = c("vector", "list"), ...)
{
	d1 = ncol(object$X)
	d2 = ncol(object$S)
	d3 = ncol(object$W)
	unfixed = object$unfixed
	idx.par1 = seq_along(unfixed$beta)
	idx.par2 = seq_along(unfixed$gamma) + length(unfixed$beta)
	idx.par3 = seq_along(unfixed$zeta) + length(unfixed$beta) +  length(unfixed$gamma)

	sd.hat = sqrt(diag(vcov(object)))

	if (match.arg(type) == "vector") {
		out = sd.hat
	} else if (match.arg(type) == "list") {
		sd.beta = rep(NA, d1)
		sd.gamma = rep(NA, d2)
		sd.zeta = rep(NA, d3)
		sd.beta[unfixed$beta] = sd.hat[idx.par1]
		sd.gamma[unfixed$gamma] = sd.hat[idx.par2]
		sd.zeta[unfixed$zeta] = sd.hat[idx.par3]
		out = list(beta = sd.beta, gamma = sd.gamma, zeta = sd.zeta)
	} else {
		stop("Unrecognized type")
	}

	return(out)
}

#' @name glm.cmp, ZICMP support
#' @export
vcov.zicmpfit = function(object, ...)
{
	# Compute the covariance via Hessian from optimizer
	-solve(object$H)
}

#' @name equitest
#' @export
equitest.zicmpfit = function(object, ...)
{
	if ("equitest" %in% names(object)) {
		return(object$equitest)
	}

	y = object$y
	X = object$X
	S = object$S
	W = object$W
	init = object$init
	offset = object$offset
	fixed = object$fixed
	control = object$control
	ll = object$loglik

	# If any elements of gamma have been fixed, an "equidispersion" test no
	# longer makes sense. Unless the values were fixed at zeroes. But let's
	# avoid this this complication.
	if (length(fixed$gamma) > 0) {
		msg = c("Chi-squared test for equidispersion not run",
			"(Some elements of gamma were fixed)")
		return(msg)
	}

	n = length(y)
	d2 = ncol(S)

	# Null model is ZICMP with nu determined by the offset off.s. If off.s happens
	# to be zeros, this simplifies to a Poisson regression.
	fit0.out = fit.zicmp.reg(y, X, S, W, offset = offset,
		init = get.init(beta = object$beta, gamma = numeric(d2), zeta = object$zeta),
		fixed = get.fixed(beta = fixed$beta, gamma = seq_len(d2), zeta = fixed$zeta),
		control = control)
	ll0 = fit0.out$loglik

	X2 = -2 * (ll0 - ll)
	pvalue = pchisq(X2, df = d2, lower.tail = FALSE)
	list(teststat = X2, pvalue = pvalue, df = d2)
}

#' @name glm.cmp, ZICMP support
#' @export
deviance.zicmpfit = function(object, ...)
{
	browser()

	y = object$y
	X = object$X
	S = object$S
	W = object$W
	n = length(y)
	# d1 = ncol(X)
	# d2 = ncol(S)
	# d3 = ncol(W)
	fixed = object$fixed
	offset = object$offset
	control = object$control

	# par.hat = c(object$beta, object$gamma, object$zeta)
	# init = get.init(beta = object$beta, gamma = object$gamma, zeta = object$zeta)
	# par.init = par.hat
	ll.star = numeric(n)

	for (i in 1:n) {
		# loglik = function(par){
		#	beta = par[1:d1]
		#	gamma = par[1:d2 + d1]
		#	zeta = par[1:d3 + d1 + d2]
		#	out = fitted.zicmp.internal(X[i,], S[i,], W[i,], beta, gamma, zeta,
		#		offset$x[i], offset$s[i], offset$w[i])
		#	dzicmp(y[i], out$lambda, out$nu, out$p, log = TRUE)
		# }

		# Maximize loglik for ith obs
		# res = optim(par.init, loglik, control = list(fnscale = -1))
		# ll.star[i] = res$value

		# Maximize loglik for ith obs
		# glm.out = glm.zicmp.raw(y[i], X[i,], S[i,], W[i,],
		#	offset = get.offset(x = offset$x[i], s = offset$s[i], w = offset$w[i]),
		#	init = par.hat, fixed = fixed, control = control)
		
		glm.out = fit.zicmp.reg(y[i],
			X = X[i,,drop = FALSE],
			S = S[i,,drop = FALSE],
			W = W[i,,drop = FALSE],
			init = get.init(beta = object$beta, gamma = numeric(d2), zeta = object$zeta),
			offset = get.offset(x = offset$x[i], s = offset$s[i], w = offset$w[i]),
			fixed = get.fixed(beta = fixed$beta, gamma = seq_len(d2), zeta = fixed$zeta),
			control = control)
		ll.star[i] = glm.out$opt.res$value
	}

	# Compute exact deviances
	ll = numeric(n)
	out = fitted.zicmp.internal(X, S, W, object$beta, object$gamma,
		object$zeta, offset$x, offset$s, offset$w)
	for (i in 1:n) {
		ll[i] = dzicmp(y[i], out$lambda[i], out$nu[i], out$p[i], log = TRUE,
			control = control)
	}

	dev = -2 * (ll - ll.star)
	leverage = leverage(object)
	cmpdev = dev / sqrt(1 - leverage)
	return(cmpdev)
}

#' @name glm.cmp, ZICMP support
#' @export
residuals.zicmpfit = function(object, type = c("raw", "quantile"), ...)
{
	out = fitted.zicmp.internal(object$X, object$S, object$W, object$beta,
		object$gamma, object$zeta, object$offset$x, object$offset$s,
		object$offset$w)

	type = match.arg(type)
	if (type == "raw") {
		y.hat = ezicmp(out$lambda, out$nu, out$p, control = object$control)
		res = object$y - y.hat
	} else if (type == "quantile") {
		res = rqres.zicmp(object$y, out$lambda, out$nu, out$p,
			control = object$control)
	} else {
		stop("Unsupported residual type")
	}

	return(as.numeric(res))
}

#' @name glm.cmp, ZICMP support
#' @export
predict.zicmpfit = function(object, newdata = NULL, type = c("response", "link"), ...)
{
	if (is.null(newdata)) {
		# If newdata is NULL, reuse data for model fit
		X = object$X
		S = object$S
		W = object$W
		off.x = object$offset$x
		off.s = object$offset$s
		off.w = object$offset$w
	} else if (object$interface == "formula") {
		# If the model was fit with the formula interface, attempt to process
		# newdata as a data.frame
		newdata = as.data.frame(newdata)

		# If the response was not included in newdata, add a column with zeros
		response.name = all.vars(object$formula.lambda)[1]
		if (is.null(newdata[[response.name]])) {
			newdata[[response.name]] = 0
		}

		raw = formula2raw(object$formula.lambda, object$formula.nu,
			object$formula.p, data = newdata, ...)
		X = raw$X
		S = raw$S
		W = raw$W
		off.x = raw$offset$x
		off.s = raw$offset$s
		off.w = raw$offset$w
	} else if (object$interface == "raw") {
		# If the model was fit with the raw interface, attempt to process
		# newdata as a list
		stopifnot(class(newdata) == "COMPoissonReg.modelmatrix")
		X = newdata$X
		S = newdata$S
		W = newdata$W
		off.x = newdata$offset$x
		off.s = newdata$offset$s
		off.w = newdata$offset$w
	} else {
		stop("Don't recognize value of interface")
	}

	link = fitted.zicmp.internal(X, S, W, object$beta,
		object$gamma, object$zeta, off.x, off.s, off.w)
	switch(match.arg(type),
		response = ezicmp(link$lambda, link$nu, link$p, control = object$control),
		link = data.frame(lambda = link$lambda, nu = link$nu, p = link$p)
	)
}

#' @name glm.cmp, ZICMP support
#' @export
parametric.bootstrap.zicmpfit = function(object, reps = 1000, report.period = reps+1, ...)
{
	n = length(object$y)
	qq = length(object$beta) + length(object$gamma) + length(object$zeta)

	out = matrix(NA, reps, qq)
	colnames(out) = c(
		sprintf("X:%s", colnames(object$X)),
		sprintf("S:%s", colnames(object$S)),
		sprintf("W:%s", colnames(object$W))
	)

	fitted.out = fitted.zicmp.internal(object$X, object$S, object$W, object$beta,
		object$gamma, object$zeta, object$offset$x, object$offset$s,
		object$offset$w)
	lambda.hat = fitted.out$lambda
	nu.hat = fitted.out$nu
	p.hat = fitted.out$p

	init = get.init(beta = object$beta, gamma = object$gamma, zeta = object$zeta)

	for (r in 1:reps) {
		if (r %% report.period == 0) {
			logger("Starting bootstrap rep %d\n", r)
		}

		# Generate bootstrap samples of the full dataset using MLE
		y.boot = rzicmp(n, lambda.hat, nu.hat, p.hat)

		# Take each of the bootstrap samples and fit model to generate bootstrap
		# estimates
		tryCatch({
			fit.boot = fit.zicmp.reg(y = y.boot, X = object$X, S = object$S,
				W = object$W, init = init, offset = object$offset,
				fixed = object$fixed, control = object$control)
			out[r,] = unlist(fit.boot$theta.hat)
		},
		error = function(e) {
			# Do nothing now; emit a warning later
		})
	}

	cnt = sum(rowSums(is.na(out)) > 0)
	if (cnt > 0) {
		warning(sprintf("%d out of %d bootstrap iterations failed", cnt, reps))
	}

	return(out)
}
