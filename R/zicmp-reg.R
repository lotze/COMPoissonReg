zicmp <- function(formula.lambda, formula.nu = NULL, formula.p = NULL,
	beta.init = NULL, gamma.init = NULL, zeta.init = NULL, max = 100, ...)
{
	# Parse formula.lambda. This one should have the response.
	# TBD: This doesn't seem robust to user passing in no response
	mf <- model.frame(formula.lambda, ...)
	y <- model.response(mf)
	X <- model.matrix(formula.lambda, mf)
	d1 <- ncol(X)

	# TBD: We should make sure that these *don't* have a response
	mf <- model.frame(formula.nu, ...)
	S <- model.matrix(formula.nu, mf)
	d2 <- ncol(S)

	# TBD: We should make sure that these *don't* have a response
	mf <- model.frame(formula.p, ...)
	W <- model.matrix(formula.p, mf)
	d3 <- ncol(W)

	initial.glm <- glm(formula.lambda, family='poisson', ...)
	if (is.null(beta.init)) { beta.init <- coef(initial.glm) }
	if (is.null(gamma.init)) { gamma.init <- rep(0, d2) }
	if (is.null(zeta.init)) { zeta.init <- rep(0, d3) }

	# TBD: Maybe include this in output object instead...
	print(beta.init)
	print(gamma.init)
	print(zeta.init)

	res <- list()
	res$formula.lambda <- formula.lambda
	res$formula.nu <- formula.nu
	res$formula.p <- formula.p
	res$y <- y
	res$X <- X
	res$S <- S
	res$W <- W
	res$max <- max

	fit.out <- fit.zicmp.reg(res$y, res$X, res$S, res$W, beta.init = beta.init,
		gamma.init = gamma.init, zeta.init = zeta.init, max = res$max)

	num_pars <- length(unlist(fit.out$theta.hat))
	res$beta.glm <- coef(initial.glm)
	res$beta <- fit.out$theta.hat$beta
	res$gamma <- fit.out$theta.hat$gamma
	res$zeta <- fit.out$theta.hat$zeta
	res$V <- fit.out$V

	attr(res, "class") <- c("zicmp", attr(res, "class"))
	return(res)
}

coef.zicmp <- function(object) {
	return(unlist(object$theta.hat))
}

nu.zicmp <- function(object, ...) {
	exp(object$S %*% object$gamma)
}

sdev.zicmp <- function(object, ...) {
	sqrt(diag(object$V))
}

chisq.zicmp <- function(object, ...) {
	LRT(object$predictors, object$response, object$glm_coefficients, object$coef, object$nu, object$max)$teststat[1,1]
}

pval.zicmp <- function(object, ...) {
	LRT(object$predictors, object$response, object$glm_coefficients, object$coef, object$nu, object$max)$pvalue
}

leverage.zicmp <- function(object, ...) {
	CMPLeverage(object$predictors, object$response, object$coef, object$nu, object$max)
}

deviance.zicmp <- function(object, ...) {
	CMPDeviance(object$predictors, object$response, object$coef, object$nu, leverage.cmp(object), object$max)
}

residuals.zicmp <- function(object, ...) {
	return(object$response - predict(object, newdata=object$predictors))
}

predict.zicmp <- function(object, ...) {
	newdata = list(...)[["newdata"]]
	return(constantCMPfitsandresids(object$coef, object$nu, newdata[,object$x_names])$fit)
}

parametric_bootstrap.zicmp <- function(object, ...) {
	n = list(...)[["n"]]
	if (is.null(n)) {
		n = 1000
	}
	bootstrap_results = as.data.frame(CMPParamBoot(x=object$predictors, object$glm_coefficients, betahat=object$coef, nuhat=object$nu, n=n)$CMPresult)
	names(bootstrap_results) = c("(Intercept)",object$x_names,"nu",recursive=TRUE)
	return(bootstrap_results)
}
