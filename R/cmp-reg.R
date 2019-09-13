#' Supporting Functions for COM-Poisson Regression
#' 
#' @param object object of type \code{cmp}.
#' @param x object of type \code{cmp}.
#' @param k Penalty per parameter to be used in AIC calculation.
#' @param newdata New covariates to be used for prediction.
#' @param type Type of residual to be computed.
#' @param reps Number of bootstrap repetitions.
#' @param report_period Report progress every \code{report_period} iterations.
#' @param ... other model parameters, such as data.
#' 
#' @name glm_cmp, CMP support
NULL

#' @name glm_cmp, CMP support
#' @export
summary.cmp = function(object, ...)
{
	n = nrow(object$X)
	d1 = ncol(object$X)
	d2 = ncol(object$S)

	V = vcov(object)
	est = coef(object)
	se = sdev(object)
	z_val = est / se
	p_val = 2*pnorm(-abs(z_val))

	DF = data.frame(
		Estimate = round(est, 4),
		SE = round(se, 4),
		z_value = round(z_val, 4),
		p_value = sprintf("%0.4g", p_val)
	)
	rownames(DF) = c(sprintf("X:%s", colnames(object$X)),
		sprintf("S:%s", colnames(object$S)))

	# If X, S, or W are intercept only, compute results for non-regression parameters
	DF_lambda = NULL
	DF_nu = NULL

	is_intercept_only = function(A, eps = 1e-12) {
		prod(dim(A) == c(n,1)) & norm(A - 1, type = "F") < eps
	}

	if (is_intercept_only(object$X)) {
		est = exp(object$beta)
		J = c(exp(object$beta), rep(0, d2))
		se = sqrt(t(J) %*% V %*% J)

		DF_lambda = data.frame(
			Estimate = round(est, 4),
			SE = round(se, 4)
		)
		rownames(DF_lambda) = "lambda"
	}

	if (is_intercept_only(object$S)) {
		est = exp(object$gamma)
		J = c(rep(0, d1), exp(object$gamma))
		se = sqrt(t(J) %*% V %*% J)

		DF_nu = data.frame(
			Estimate = round(est, 4),
			SE = round(se, 4)
		)
		rownames(DF_nu) = "nu"
	}

	list(DF = DF, DF_lambda = DF_lambda, DF_nu = DF_nu,
		n = n,
		loglik = logLik(object),
		aic = AIC(object),
		bic = BIC(object),
		opt_method = object$opt_method,
		opt_message = object$opt_res$message,
		opt_convergence = object$opt_res$convergence,
		elapsed_sec = object$elapsed_sec
	)
}

fitted_cmp_internal = function(X, S, beta, gamma, off_x, off_s)
{
	list(
		lambda = as.numeric(exp(X %*% beta + off_x)),
		nu = as.numeric(exp(S %*% gamma + off_s))
	)
}

#' @name glm_cmp, CMP support
#' @export
print.cmp = function(x, ...)
{
	printf("CMP coefficients\n")
	s = summary.cmp(x)
	tt = equitest(x)
	print(s$DF)

	if (!is.null(s$DF_lambda) || !is.null(s$DF_nu)) {
		printf("--\n")
		printf("Transformed intercept-only parameters (excluding offsets)\n")
		print(rbind(s$DF_lambda, s$DF_nu))
	}
	printf("--\n")
	printf("Chi-squared test for equidispersion\n")
	printf("X^2 = %0.4f, df = 1, ", tt$teststat)
	printf("p-value = %0.4e\n", tt$pvalue)
	printf("--\n")
	printf("Elapsed Sec: %0.2f   ", s$elapsed_sec)
	printf("Sample size: %d   ", s$n)
	printf("SEs via Hessian\n")
	printf("LogLik: %0.4f   ", s$loglik)
	printf("AIC: %0.4f   ", s$aic)
	printf("BIC: %0.4f   ", s$bic)
	printf("\n")
	printf("Optimization Method: %s   ", s$opt_method)
	printf("Converged status: %d\n", s$opt_convergence)
	printf("Message: %s\n", s$opt_message)
}

#' @name glm_cmp, CMP support
#' @export
logLik.cmp = function(object, ...)
{
	object$loglik
}

#' @name glm_cmp, CMP support
#' @export
AIC.cmp = function(object, ..., k=2)
{
	-2*object$loglik + 2*length(coef(object))
}

#' @name glm_cmp, CMP support
#' @export
BIC.cmp = function(object, ...)
{
	n = length(object$y)
	-2*object$loglik + log(n)*length(coef(object))
}

#' @name glm_cmp, CMP support
#' @export
coef.cmp = function(object, ...)
{
	c(object$beta, object$gamma)
}

#' @name glm_cmp, CMP support
#' @export
nu.cmp = function(object, ...)
{
	out = fitted_cmp_internal(object$X, object$S, object$beta, object$gamma,
		object$off_x, object$off_s)
	out$nu
}

#' @name glm_cmp, CMP support
#' @export
sdev.cmp = function(object, ...)
{
	sqrt(diag(vcov(object)))
}

#' @name glm_cmp, CMP support
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

	out = fitted_cmp_internal(object$X, object$S, object$beta, object$gamma,
		object$off_x, object$off_s)
	lambda_hat = out$lambda
	nu_hat = out$nu

	out0 = fitted_cmp_internal(object$X, object$S, object$beta_glm,
		gamma = numeric(d2), object$off_x, object$off_s)
	lambda0_hat = out0$lambda

	logz = z_hybrid(lambda_hat, nu_hat, take_log = TRUE)
	teststat = -2 * sum(y*log(lambda0_hat) - lgamma(y+1) - lambda0_hat -
		y*log(lambda_hat) + nu_hat*lgamma(y+1) + logz)
	pvalue = pchisq(teststat, df = 1, lower.tail = FALSE)
	list(teststat = teststat, pvalue = pvalue)
}

#' @name glm_cmp, CMP support
#' @export
leverage.cmp = function(object, ...)
{
	y = object$y

	# 1) to code the WW matrix  (diagonal matrix with Var(Y_i) )
	ww = weights(object$X, object$S, object$beta, object$gamma, object$off_x, object$off_s)
	WW = diag(ww)

	out = fitted_cmp_internal(object$X, object$S, object$beta, object$gamma,
		object$off_x, object$off_s)
	lambda_hat = out$lambda
	nu_hat = out$nu

	#    and X matrix (in Appendix)
	E_y = z_prodj(lambda_hat, nu_hat) / z_hybrid(lambda_hat, nu_hat)
	E_logfacty = z_prodlogj(lambda_hat, nu_hat) / z_hybrid(lambda_hat, nu_hat)
	extravec = (-lgamma(y+1) + E_logfacty)/(y - E_y)
	curlyX_mat = cbind(object$X, extravec)

	# 2) to compute H using eq (12)  on p. 11
	H1 = t(curlyX_mat) %*% sqrt(WW)
	H2 = solve(t(curlyX_mat) %*% WW %*% curlyX_mat)
	H = t(H1) %*% H2 %*% H1
	diagH = diag(H)
	return(diagH)
}

#' @name glm_cmp, CMP support
#' @export
deviance.cmp = function(object, ...)
{
	# Compute the COM-Poisson deviances exactly
	y = object$y
	n = length(y)

	#### Compute optimal log likelihood value for given nu-hat value
	beta_init = object$beta
	d1 = length(beta_init)
	ll_y = numeric(n)

	for (i in 1:n) {
		# loglik for single observation
		logf = function(beta) {
			out = fitted_cmp_internal(object$X[i,], object$S[i,], beta, object$gamma,
				object$off_x[i], object$off_s[i])
			y[i]*log(out$lambda) - out$nu*lgamma(y[i] + 1) -
				z_hybrid(out$lambda, out$nu, take_log = TRUE)
		}

		# Determine the MLEs
		# Using optim rather than nlm because optim handles -Inf more gracefully, if encountered
		opt_method = getOption("COMPoissonReg.optim.method")
		res = optim(beta_init, logf, method = opt_method, control = list(fnscale = -1))
		ll_y[i] = res$value
	}

	# Compute exact deviances
	out = fitted_cmp_internal(object$X, object$S, object$beta, object$gamma,
		object$off_x, object$off_s)
	ll_mu = y*log(out$lambda) - out$nu * lgamma(y+1) -
		z_hybrid(out$lambda, out$nu, take_log = TRUE)
	d = -2*(ll_mu - ll_y)
	lev = leverage.cmp(object)
	cmpdev = d / sqrt(1 - lev)
	return(cmpdev)
}

#' @name glm_cmp, CMP support
#' @export
residuals.cmp = function(object, type = c("raw", "quantile"), ...)
{
	out = fitted_cmp_internal(object$X, object$S, object$beta, object$gamma, object$off_x, object$off_s)
	y_hat = predict.cmp(object, newdata = object$X)

	type = match.arg(type)
	if (type == "raw") {
		res = object$y - y_hat
	} else if (type == "quantile") {
		res = rqres.cmp(object$y, lambda = out$lambda, nu = out$nu)
	} else {
		stop("Unsupported residual type")
	}

	return(as.numeric(res))
}

#' @name glm_cmp, CMP support
#' @export
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

	out = fitted_cmp_internal(X, S, object$beta, object$gamma, object$off_x, object$off_s)
	fnr = constantCMPfitsandresids(out$lambda, out$nu)
	y_hat = fnr$fit
	return(y_hat)
}

#' @name glm_cmp, CMP support
#' @export
parametric_bootstrap.cmp = function(object, reps = 1000, report_period = reps+1, ...)
{
	n = nrow(object$X)
	d1 = ncol(object$X)
	d2 = ncol(object$S)

	out = fitted_cmp_internal(object$X, object$S, object$beta, object$gamma, object$off_x, object$off_s)
	lambda_hat = out$lambda
	nu_hat = out$nu

	# Generate `reps` samples, using beta_hat and nu_hat from full dataset
	# Run CMP regression on each boostrap sample to generate new beta and nu estimates
	
	boot_out = matrix(NA, nrow = reps, ncol = d1 + d2)

	for (r in 1:reps){
		if (r %% report_period == 0) {
			logger("Starting bootstrap rep %d\n", r)
		}
		y_boot = rcmp(n, lambda_hat, nu_hat)
		tryCatch({
			res = fit_cmp_reg(y = y_boot, X = object$X, S = object$S,
				beta_init = object$beta_glm, gamma_init = object$gamma,
				off_x = object$off_x, off_s = object$off_s)
			boot_out[r,] = unlist(res$theta_hat)
		}, error = function(e) {
			# Do nothing now; emit a warning later
		})
	}

	cnt = sum(rowSums(is.na(boot_out)) > 0)
	if (cnt > 0) {
		warning(sprintf("%d out of %d bootstrap iterations failed", cnt, reps))
	}

	colnames(boot_out) = c(colnames(object$X), colnames(object$S), recursive=TRUE)
	return(boot_out)
}

constantCMPfitsandresids = function(lambda_hat, nu_hat, y = 0)
{
	# Determine estimated lambda_hat, fit, and residuals
	lambda_hat = as.numeric(lambda_hat)
	nu_hat = as.numeric(nu_hat)
	fit = lambda_hat^(1/nu_hat) - ((nu_hat - 1)/(2*nu_hat))
	resid = y - fit
	list(fit = fit, resid = resid)
}

weights = function(X, S, beta, gamma, off_x, off_s)
{
	out = fitted_cmp_internal(X, S, beta, gamma, off_x, off_s)

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

cmp_mse = function(y, X, S, beta, gamma, off_x, off_s)
{
	out = fitted_cmp_internal(X, S, beta, gamma, off_x, off_s)
	fnr = constantCMPfitsandresids(out$lambda, out$nu, y)
	res = fnr$resid
	mse = mean(res^2)
	return(mse)
}
