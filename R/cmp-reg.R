#' Supporting Functions for COM-Poisson Regression
#' 
#' @param object object of type \code{cmp}.
#' @param x object of type \code{cmp}.
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
#' columns corresponding to estimates of \code{lambda} and \code{nu}.
#' 
#' The function \code{coef} returns a vector of coefficient estimates in
#' the form \code{c(beta, gamma)} when \code{type = "vector"}. When
#' \code{type = "list"}, the estimates are returned as a list with named
#' elements \code{beta} and \code{gamma}.
#' 
#' The \code{type} argument behaves the same for the \code{sdev} function
#' as it does for \code{coef}.
#' 
#' @name glm.cmp, CMP support
NULL

#' @name glm.cmp, CMP support
#' @export
summary.cmpfit = function(object, ...)
{
	n = nrow(object$X)
	d1 = ncol(object$X)
	d2 = ncol(object$S)
	qq = d1 + d2

	# We need the indices of the fixed coefficients and the ones included in
	# optimization.
	fixed = object$fixed
	unfixed = object$unfixed
	idx.par1 = seq_along(unfixed$beta)
	idx.par2 = seq_along(unfixed$gamma) + length(unfixed$beta)

	V = vcov(object)
	est = coef(object)

	# In the vector of SEs, include an NA entry if the variable was fixed
	# This NA will propagate to the corresponding z-value and p-value as well.
	se = rep(NA, qq)
	se[unfixed$beta] = sdev(object)[idx.par1]
	se[unfixed$gamma + d1] = sdev(object)[idx.par2]

	z.val = est / se
	p.val = 2*pnorm(-abs(z.val))

	DF = data.frame(
		Estimate = round(est, 4),
		SE = round(se, 4),
		z.value = round(z.val, 4),
		p.value = sprintf("%0.4g", p.val)
	)
	colnames(DF) = c("Estimate", "SE", "z-value", "p-value")
	rownames(DF) = c(
		sprintf("X:%s", colnames(object$X)),
		sprintf("S:%s", colnames(object$S))
	)

	# Add a column to indicate fixed components if anything is fixed
	if (length(unlist(fixed)) > 0) {
		fixed.beta = rep("F", d1)
		fixed.gamma = rep("F", d2)
		fixed.beta[fixed$beta] = "T"
		fixed.gamma[fixed$gamma] = "T"
		DF$Fixed = c(fixed.beta, fixed.gamma)
	}

	# If X, S, or W are intercept only, compute results for non-regression parameters
	DF.lambda = NULL
	DF.nu = NULL

	# In each block below, make sure to consider only the non-fixed variables
	# for the Jacobian and Hessian. If one of the intercepts was fixed, it should
	# result in an SE of zero.

	if (is.intercept.only(object$X) && is.zero.matrix(object$offset$x)) {
		if (length(fixed$beta) > 0) {
			se = 0
		} else {
			J = c(exp(object$beta), numeric(length(idx.par2)))
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
			J = c(numeric(length(idx.par1)), exp(object$gamma))
			se = sqrt(t(J) %*% V %*% J)
		}
		est = exp(object$gamma)
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
		optim.method = object$control$optim.method,
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
print.cmpfit = function(x, ...)
{
	printf("CMP coefficients\n")
	s = summary(x)
	print(s$DF)
	
	if (length(x$fixed$gamma) > 0) {
		tt = paste(collapse = " ", c(
			"Some elements of gamma were fixed.",
			"Chi-squared test for equidispersion not defined."))
	} else {
		tt = equitest(x)
	}

	if (!is.null(s$DF.lambda) || !is.null(s$DF.nu)) {
		printf("--\n")
		printf("Transformed intercept-only parameters\n")
		print(rbind(s$DF.lambda, s$DF.nu))
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

#' @name glm.cmp, CMP support
#' @export
logLik.cmpfit = function(object, ...)
{
	object$loglik
}

#' @name glm.cmp, CMP support
#' @export
AIC.cmpfit= function(object, ..., k=2)
{
	d = length(unlist(object$unfixed))
	-2*object$loglik + 2*d
}

#' @name glm.cmp, CMP support
#' @export
BIC.cmpfit = function(object, ...)
{
	n = length(object$y)
	d = length(unlist(object$unfixed))
	-2*object$loglik + log(n)*d
}

#' @name glm.cmp, CMP support
#' @export
coef.cmpfit = function(object, type = c("vector", "list"), ...)
{
	switch(match.arg(type),
		vector = c(object$beta, object$gamma),
		list = list(beta = object$beta, gamma = object$gamma)
	)
}

#' @name glm.cmp, CMP support
#' @export
nu.cmpfit = function(object, ...)
{
	# This function is deprecated - use predict instead
	.Deprecated("predict(object, type = \"link\")")
	link = predict(object, type = "link")
	link$nu
}

#' @name glm.cmp, CMP support
#' @export
sdev.cmpfit = function(object, type = c("vector", "list"), ...)
{
	d1 = ncol(object$X)
	d2 = ncol(object$S)
	fixed = object$fixed
	unfixed = object$unfixed
	idx.par1 = seq_along(unfixed$beta)
	idx.par2 = seq_along(unfixed$gamma) + length(unfixed$beta)

	sd.hat = sqrt(diag(vcov(object)))

	if (match.arg(type) == "vector") {
		out = sd.hat
	} else if (match.arg(type) == "list") {
		sd.beta = rep(NA, d1)
		sd.gamma = rep(NA, d2)
		sd.beta[unfixed$beta] = sd.hat[idx.par1]
		sd.gamma[unfixed$gamma] = sd.hat[idx.par2]
		out = list(beta = sd.beta, gamma = sd.gamma)
	} else {
		stop("Unrecognized type")
	}

	return(out)
}

#' @name glm.cmp, CMP support
#' @export
vcov.cmpfit = function(object, ...)
{
	# Compute the covariance via Hessian from optimizer
	-solve(object$H)
}

#' @name equitest
#' @export
equitest.cmpfit = function(object, ...)
{
	if ("equitest" %in% names(object)) {
		return(object$equitest)
	}

	y = object$y
	X = object$X
	S = object$S
	init = object$init
	offset = object$offset
	fixed = object$fixed
	ll = object$loglik

	# If any elements of gamma have been fixed, an "equidispersion" test no
	# longer makes sense. Unless the values were fixed at zeroes. But let's
	# avoid this complication.
	if (length(fixed$gamma) > 0) {
		msg = c("Some elements of gamma were fixed,",
			"chi-squared test for equidispersion not defined")
		stop(paste(msg, collapse = " "))
	}

	n = length(y)
	d2 = ncol(S)

	# Null model is CMP with nu determined by the offset off.s. If off.s happens
	# to be zeros, this simplifies to a Poisson regression.
	fit0.out = fit.cmp.reg(y, X, S, offset = offset,
		init = get.init(beta = object$beta, gamma = numeric(d2), zeta = numeric(0)),
		fixed = get.fixed(beta = fixed$beta, gamma = seq_len(d2), zeta = fixed$zeta),
		control = object$control)
	ll0 = fit0.out$loglik

	teststat = -2 * (ll0 - ll)
	pvalue = pchisq(teststat, df = d2, lower.tail = FALSE)
	list(teststat = teststat, pvalue = pvalue, df = d2)
}

#' @name glm.cmp, CMP support
#' @export
leverage.cmpfit = function(object, ...)
{
	y = object$y
	n = length(y)

	out = fitted.cmp.internal(object$X, object$S, object$beta, object$gamma,
		object$offset$x, object$offset$s)

	# 1) Some quantities corresponding to parameters lambda and nu: Normalizing
	# constants, expected values, variances, and truncation values.
	z.hat = ncmp(out$lambda, out$nu, control = object$control)
	E.y = ecmp(out$lambda, out$nu, control = object$control)
	V.y = vcmp(out$lambda, out$nu, control = object$control)
	y.trunc = max(tcmp(out$lambda, out$nu, control = object$control))

	# Note that z_prodlogj uses truncation, even if we would approximate z 
	# using the asymptotic expression using some of the given lambda and nu
	# pairs. This might be something that could be improved later.

	#    and X matrix (in Appendix)
	E.logfacty = z_prodlogj(out$lambda, out$nu, max = y.trunc) / z.hat
	extravec = (E.logfacty - lgamma(y+1)) / (y - E.y)
	curlyX.mat = cbind(object$X, extravec)

	# 2) to compute H using equation (12)  on p. 11
	if (FALSE) {
		# Here is a more readable version, but it creates n x n matrices so may
		# have problems with large inputs.
		WW = diag(V.y)
		H1 = t(curlyX.mat) %*% sqrt(WW)
		H2 = solve(t(curlyX.mat) %*% WW %*% curlyX.mat)
		H = t(H1) %*% H2 %*% H1
		diag.H = diag(H)
	} else {
		# This version avoids creating large matrices, but is a bit harder to read.
		diag.H = numeric(n)
		H1 = t(curlyX.mat * sqrt(V.y))
		H2 = solve(t(curlyX.mat * V.y) %*% curlyX.mat)
		for (i in 1:n) {
			diag.H[i] = t(H1[,i]) %*% H2 %*% H1[,i]
		}
	}

	return(diag.H)
}

#' @name glm.cmp, CMP support
#' @export
deviance.cmpfit = function(object, ...)
{
	# Compute the COM-Poisson deviances exactly
	y = object$y
	X = object$X
	S = object$S
	# init = object$init
	fixed = object$fixed
	offset = object$offset
	control = object$control
	n = length(y)
	d2 = ncol(S)

	# Compute optimal log likelihood value for given nu-hat value
	# beta.init = object$beta
	# d1 = length(beta.init)
	ll.star = numeric(n)

	for (i in 1:n) {
		# Maximize loglik for ith obs
		glm.out = fit.cmp.reg(y[i],
			X = X[i,,drop = FALSE],
			S = S[i,,drop = FALSE],
			init = get.init(beta = object$beta, gamma = numeric(d2)),
			offset = get.offset(x = offset$x[i], s = offset$s[i], w = offset$w[i]),
			fixed = get.fixed(beta = fixed$beta, gamma = seq_len(d2), zeta = fixed$zeta),
			control = control)
		ll.star[i] = glm.out$opt.res$value
	}

	# Compute exact deviances
	ll = numeric(n)
	out = fitted.cmp.internal(X, S, object$beta, object$gamma, offset$x, offset$s)
	for (i in 1:n) {
		ll[i] = dcmp(y[i], out$lambda[i], out$nu[i], log = TRUE, control = control)
	}

	dev = -2*(ll - ll.star)
	lev = leverage(object)
	cmpdev = dev / sqrt(1 - lev)
	return(cmpdev)
}

#' @name glm.cmp, CMP support
#' @export
residuals.cmpfit = function(object, type = c("raw", "quantile"), ...)
{
	out = fitted.cmp.internal(object$X, object$S, object$beta, object$gamma,
		object$offset$x, object$offset$s)

	type = match.arg(type)
	if (type == "raw") {
		y.hat = ecmp(out$lambda, out$nu, control = object$control)
		res = object$y - y.hat
	} else if (type == "quantile") {
		res = rqres.cmp(object$y, lambda = out$lambda, nu = out$nu,
			control = object$control)
	} else {
		stop("Unsupported residual type")
	}

	return(as.numeric(res))
}

#' @name glm.cmp, CMP support
#' @export
predict.cmpfit = function(object, newdata = NULL, type = c("response", "link"), ...)
{
	if (is.null(newdata)) {
		# If newdata is NULL, reuse data from model fit
		X = object$X
		S = object$S
		off.x = object$offset$x
		off.s = object$offset$s
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
		off.x = raw$offset$x
		off.s = raw$offset$s
	} else if (object$interface == "raw") {
		# If the model was fit with the raw interface, attempt to process
		# newdata as a list
		stopifnot(class(newdata) == "COMPoissonReg.modelmatrix")
		X = newdata$X
		S = newdata$S
		off.x = newdata$offset$x
		off.s = newdata$offset$s
	} else {
		stop("Don't recognize value of interface")
	}

	link = fitted.cmp.internal(X, S, object$beta, object$gamma, off.x, off.s)
	switch(match.arg(type),
		response = ecmp(link$lambda, link$nu),
		link = data.frame(lambda = link$lambda, nu = link$nu)
	)
}

#' @name glm.cmp, CMP support
#' @export
parametric.bootstrap.cmpfit = function(object, reps = 1000, report.period = reps+1, ...)
{
	n = nrow(object$X)
	d1 = ncol(object$X)
	d2 = ncol(object$S)
	qq = d1 + d2

	out = matrix(NA, nrow = reps, ncol = qq)
	colnames(out) = c(colnames(object$X), colnames(object$S), recursive=TRUE)

	fitted.out = fitted.cmp.internal(object$X, object$S, object$beta, object$gamma,
		object$offset$x, object$offset$s)
	lambda.hat = fitted.out$lambda
	nu.hat = fitted.out$nu

	init = get.init(beta = object$beta, gamma = object$gamma)

	for (r in 1:reps){
		if (r %% report.period == 0) {
			logger("Starting bootstrap rep %d\n", r)
		}

		# Generate bootstrap samples of the full dataset using MLE
		y.boot = rcmp(n, lambda.hat, nu.hat, control = object$control)

		# Take each of the bootstrap samples and fit model to generate bootstrap
		# estimates
		tryCatch({
			fit.boot = fit.cmp.reg(y = y.boot, X = object$X, S = object$S,
				init = init, offset = object$offset, fixed = object$fixed,
				control = object$control)
			out[r,] = unlist(fit.boot$theta.hat)
		}, error = function(e) {
			# Do nothing now; emit a warning later
		})
	}

	cnt = sum(rowSums(is.na(out)) > 0)
	if (cnt > 0) {
		warning(sprintf("%d out of %d bootstrap iterations failed", cnt, reps))
	}

	return(out)
}
