#' Supporting Functions for ZICMP Regression
#' 
#' @param object object of type \code{zicmp}.
#' @param x object of type \code{zicmp}.
#' @param k Penalty per parameter to be used in AIC calculation.
#' @param newdata New covariates to be used for prediction.
#' @param type Type of residual to be computed.
#' @param reps Number of bootstrap repetitions.
#' @param report_period Report progress every \code{report_period} iterations.
#' @param ... other model parameters, such as data.
#' 
#' @name glm_cmp, ZICMP support
NULL

#' @name glm_cmp, ZICMP support
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
	z_val = est / se
	p_val = 2*pnorm(-abs(z_val))
	qq = length(est)

	DF = data.frame(
		Estimate = round(est, 4),
		SE = round(se, 4),
		z_value = round(z_val, 4),
		p_value = sprintf("%0.4g", p_val)
	)
	rownames(DF) = c(
		sprintf("X:%s", colnames(object$X)),
		sprintf("S:%s", colnames(object$S)),
		sprintf("W:%s", colnames(object$W))
	)

	# If X, S, or W are intercept only, compute results for non-regression parameters
	# Exclude offsets from these calculations
	DF_lambda = NULL
	DF_nu = NULL
	DF_p = NULL

	is_intercept_only = function(A, eps = 1e-12) {
		prod(dim(A) == c(n,1)) & norm(A - 1, type = "F") < eps
	}

	if (is_intercept_only(object$X)) {
		est = exp(object$beta)
		J = c(exp(object$beta), rep(0, d2), rep(0, d3))
		se = sqrt(t(J) %*% V %*% J)

		DF_lambda = data.frame(
			Estimate = round(est, 4),
			SE = round(se, 4)
		)
		rownames(DF_lambda) = "lambda"
	}

	if (is_intercept_only(object$S)) {
		est = exp(object$gamma)
		J = c(rep(0, d1), exp(object$gamma), rep(0, d3))
		se = sqrt(t(J) %*% V %*% J)

		DF_nu = data.frame(
			Estimate = round(est, 4),
			SE = round(se, 4)		)
		rownames(DF_nu) = "nu"
	}

	if (is_intercept_only(object$W)) {
		est = plogis(object$zeta)
		J = c(rep(0, d1), rep(0, d2), dlogis(object$zeta))
		se = sqrt(t(J) %*% V %*% J)

		DF_p = data.frame(
			Estimate = round(est, 4),
			SE = round(se, 4)
		)
		rownames(DF_p) = "p"
	}

	list(DF = DF, DF_lambda = DF_lambda, DF_nu = DF_nu, DF_p = DF_p,
		n = length(object$y),
		loglik = logLik(object),
		aic = AIC(object),
		bic = BIC(object),
		opt_method = object$opt_method,
		opt_res = object$opt_res,
		elapsed_sec = object$elapsed_sec
	)
}

fitted_zicmp_internal = function(X, S, W, beta, gamma, zeta, off_x, off_s, off_w)
{
	list(
		lambda = as.numeric(exp(X %*% beta + off_x)),
		nu = as.numeric(exp(S %*% gamma + off_s)),
		p = as.numeric(plogis(W %*% zeta + off_w))
	)
}

#' @name glm_cmp, ZICMP support
#' @export
print.zicmp = function(x, ...)
{
	printf("ZICMP coefficients\n")
	s = summary(x)
	tt = equitest(x)
	print(s$DF)

	if (!is.null(s$DF_lambda) || !is.null(s$DF_nu) || !is.null(s$DF_p)) {
		printf("--\n")
		printf("Transformed intercept-only parameters (excluding offsets)\n")
		print(rbind(s$DF_lambda, s$DF_nu, s$DF_p))
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
	printf("Converged status: %d   ", s$opt_convergence)
	printf("Message: %s\n", s$opt_message)
}

#' @name glm_cmp, ZICMP support
#' @export
logLik.zicmp = function(object, ...)
{
	object$loglik
}

#' @name glm_cmp, ZICMP support
#' @export
AIC.zicmp = function(object, ..., k = 2)
{
	-2*object$loglik + 2*length(coef(object))
}

#' @name glm_cmp, ZICMP support
#' @export
BIC.zicmp = function(object, ...)
{
	n = length(object$y)
	-2*object$loglik + log(n)*length(coef(object))
}

#' @name glm_cmp, ZICMP support
#' @export
coef.zicmp = function(object, ...)
{
	c(object$beta, object$gamma, object$zeta)
}

#' @name glm_cmp, ZICMP support
#' @export
nu.zicmp = function(object, ...)
{
	out = fitted_zicmp_internal(object$X, object$S, object$W, object$beta,
		object$gamma, object$zeta, object$off_x, object$off_s, object$off_w)
	out$nu
}

#' @name glm_cmp, ZICMP support
#' @export
sdev.zicmp = function(object, ...)
{
	sqrt(diag(vcov(object)))
}

#' @name glm_cmp, ZICMP support
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
	beta_hat = object$beta
	gamma_hat = object$gamma
	zeta_hat = object$zeta
	off_x = object$off_x
	off_s = object$off_s
	off_w = object$off_w
	n = length(y)
	d1 = ncol(X)
	d2 = ncol(S)
	d3 = ncol(W)

	fit0.out = fit_zip_reg(y, X, W, beta_init = beta_hat,
		zeta_init = zeta_hat, off_x = off_x, off_w = off_w)
	beta0_hat = fit0.out$theta_hat$beta
	zeta0_hat = fit0.out$theta_hat$zeta

	out1 = fitted_zicmp_internal(X, S, W, beta_hat, gamma_hat, zeta_hat, off_x, off_s, off_w)
	lambda_hat = out1$lambda
	nu_hat = out1$nu
	p_hat = out1$p

	out0 = fitted_zicmp_internal(X, S, W, beta0_hat, numeric(d2), zeta0_hat, off_x, off_s, off_w)
	lambda0_hat = out0$lambda
	p0_hat = out0$p

	logff = dzicmp(y, lambda_hat, nu_hat, p_hat, log = TRUE)
	logff0 = dzip(y, lambda0_hat, p0_hat, log = TRUE)

	X2 = 2*(sum(logff) - sum(logff0))
	pvalue = pchisq(X2, df = d3, lower.tail = FALSE)
	list(teststat = X2, pvalue = pvalue)
}

#' @name glm_cmp, ZICMP support
#' @export
leverage.zicmp = function(object, ...)
{
	stop("This function is not yet implemented")
}

#' @name glm_cmp, ZICMP support
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
	off_x = object$off_x
	off_s = object$off_s
	off_w = object$off_w

	par_hat = c(object$beta, object$gamma, object$zeta)
	par_init = par_hat
	ll = numeric(n)
	ll_star = numeric(n)

	for (i in 1:n) {
		loglik = function(par){
			beta = par[1:d1]
			gamma = par[1:d2 + d1]
			zeta = par[1:d3 + d1 + d2]
			out = fitted_zicmp_internal(X[i,], S[i,], W[i,], beta, gamma, zeta,
				off_x[i], off_s[i], off_w[i])
			dzicmp(y[i], out$lambda, out$nu, out$p, log = TRUE)
		}

		# Maximize loglik for ith obs
		res = optim(par_init, loglik, control = list(fnscale = -1))
		ll_star[i] = res$value

		# loglik maximized over all the obs, evaluated for ith obs
		ll[i] = loglik(par_hat)
	}

	dev = -2*(ll - ll_star)
	leverage = leverage(object)
	cmpdev = dev / sqrt(1 - leverage)
	return(cmpdev)
}

#' @name glm_cmp, ZICMP support
#' @export
residuals.zicmp = function(object, type = c("raw", "quantile"), ...)
{
	out = fitted_zicmp_internal(object$X, object$S, object$W, object$beta,
		object$gamma, object$zeta, object$off_x, object$off_s, object$off_w)
	lambda_hat = out$lambda
	nu_hat = out$nu
	p_hat = out$p
	y_hat = predict.zicmp(object)

	type = match.arg(type)
	if (type == "raw") {
		res = object$y - y_hat
	} else if (type == "quantile") {
		res = rqres.zicmp(object$y, lambda_hat, nu_hat, p_hat)
	} else {
		stop("Unsupported residual type")
	}

	return(as.numeric(res))
}

#' @name glm_cmp, ZICMP support
#' @export
predict.zicmp = function(object, newdata = NULL, ...)
{
	if (!is.null(newdata)) {
		# If any of the original models had an intercept added via model.matrix, they
		# will have an "(Intercept)" column. Let's add an "(Intercept)" to newdata
		# in case the user didn't make one.
		newdata$'(Intercept)' = 1

		X = as.matrix(newdata[,colnames(object$X)])
		S = as.matrix(newdata[,colnames(object$S)])
		W = as.matrix(newdata[,colnames(object$W)])
	} else {
		X = object$X
		S = object$S
		W = object$W
	}

	out = fitted_zicmp_internal(X, S, W, object$beta,
		object$gamma, object$zeta, object$off_x, object$off_s, object$off_w)
	y_hat = zicmp_expected_value(out$lambda, out$nu, out$p)
	return(y_hat)
}

#' @name glm_cmp, ZICMP support
#' @export
parametric_bootstrap.zicmp = function(object, reps = 1000, report_period = reps+1, ...)
{
	n = length(object$y)
	qq = length(object$beta) + length(object$gamma) + length(object$zeta)
	theta_boot = matrix(NA, reps, qq)

	out = fitted_zicmp_internal(object$X, object$S, object$W, object$beta,
		object$gamma, object$zeta, object$off_x, object$off_s, object$off_w)
	lambda_hat = out$lambda
	nu_hat = out$nu
	p_hat = out$p

	for (r in 1:reps) {
		if (r %% report_period == 0) {
			logger("Starting bootstrap rep %d\n", r)
		}

		# Generate bootstrap samples of the full dataset using MLE
		y_boot = rzicmp(n, lambda_hat, nu_hat, p_hat)

		# Take each of the bootstrap samples and fit model to generate bootstrap
		# estimates
		tryCatch({
			fit_boot = fit_zicmp_reg(y_boot, object$X, object$S, object$W,
				object$beta_init, object$gamma_init, object$zeta_init,
				off_x = object$off_x, off_s = object$off_s, off_w = object$off_w)
			theta_boot[r,] = unlist(fit_boot$theta_hat)
		},
		error = function(e) {
			# Do nothing now; emit a warning later
		})
	}

	cnt = sum(rowSums(is.na(theta_boot)) > 0)
	if (cnt > 0) {
		warning(sprintf("%d out of %d bootstrap iterations failed", cnt, reps))
	}

	colnames(theta_boot) = c(
		sprintf("X:%s", colnames(object$X)),
		sprintf("S:%s", colnames(object$S)),
		sprintf("W:%s", colnames(object$W))
	)

	return(theta_boot)
}
