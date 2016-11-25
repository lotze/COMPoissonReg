zicmp <- function(formula.lambda, formula.nu = ~ 1, formula.p = NULL,
	beta.init = NULL, gamma.init = NULL, zeta.init = NULL, max = 100, ...)
{
	# Parse formula.lambda. This one should have the response.
	mf <- model.frame(formula.lambda, ...)
	y <- model.response(mf)
	X <- model.matrix(formula.lambda, mf)
	d1 <- ncol(X)

	# Parse formula.nu
	mf <- model.frame(formula.nu, ...)
	S <- model.matrix(formula.nu, mf)
	d2 <- ncol(S)

	# Parse formula.p
	mf <- model.frame(formula.p, ...)
	W <- model.matrix(formula.p, mf)
	d3 <- ncol(W)

	initial.glm <- glm(formula.lambda, family='poisson', ...)
	if (is.null(beta.init)) { beta.init <- coef(initial.glm) }
	if (is.null(gamma.init)) { gamma.init <- rep(0, d2) }
	if (is.null(zeta.init)) { zeta.init <- rep(0, d3) }

	res <- list()
	res$formula.lambda <- formula.lambda
	res$formula.nu <- formula.nu
	res$formula.p <- formula.p
	res$y <- y
	res$X <- X
	res$S <- S
	res$W <- W
	res$max <- max

	# TBD: Maybe include this in output object instead...
	res$beta.init <- beta.init
	res$gamma.init <- gamma.init
	res$zeta.init <- zeta.init

	fit.out <- fit.zicmp.reg(res$y, res$X, res$S, res$W, beta.init = beta.init,
		gamma.init = gamma.init, zeta.init = zeta.init, max = res$max)

	res$beta.glm <- coef(initial.glm)
	res$beta <- fit.out$theta.hat$beta
	res$gamma <- fit.out$theta.hat$gamma
	res$zeta <- fit.out$theta.hat$zeta
	res$V <- fit.out$V
	res$loglik <- fit.out$loglik
	res$opt.res <- fit.out$opt.res
	res$elapsed.sec <- fit.out$elapsed.sec

	attr(res, "class") <- c("zicmp", attr(res, "class"))
	return(res)
}

summary.zicmp <- function(object, ...)
{
	est <- coef(object)
	se <- sdev(object)
	z.val <- est / se
	p.val <- 2*(1 - pnorm(abs(est / se)))

	DF <- data.frame(
		Estimate = round(est, 4),
		SE = round(se, 4),
		z.value = round(z.val, 6),
		p.value = sprintf("%0.4g", p.val)
	)
	rownames(DF) <- c(
		sprintf("X:%s", colnames(object$X)),
		sprintf("S:%s", colnames(object$S)),
		sprintf("W:%s", colnames(object$W))
	)

	list(DF = DF,
		n = length(object$y),
		loglik = logLik(object),
		aic = AIC(object),
		bic = BIC(object),
		opt.res = object$opt.res,
		elapsed.sec = object$elapsed.sec
	)
}

print.zicmp <- function(x, ...)
{
	cat("Fit for ZICMP model\n")
	s <- summary(x)
	print(s$DF)

	cat("--\n")
	cat(sprintf("Elapsed Sec: %0.2f   ", s$elapsed.sec))
	cat(sprintf("Sample size: %d\n", s$n))
	cat(sprintf("LogLik: %0.4f   ", s$loglik))
	cat(sprintf("AIC: %0.4f   ", s$aic))
	cat(sprintf("BIC: %0.4f   ", s$bic))
	cat("\n")
	cat(sprintf("Converged status: %d   ", s$opt.res$convergence))
	cat(sprintf("Message: %s\n", s$opt.res$message))
}

logLik.zicmp <- function(object, k, ...)
{
	object$loglik
}

AIC.zicmp <- function(object, k, ...)
{
	-2*object$loglik + 2*length(coef(object))
}

BIC.zicmp <- function(object, ...)
{
	n <- length(object$y)
	-2*object$loglik + log(n)*length(coef(object))
}

coef.zicmp <- function(object)
{
	c(object$beta, object$gamma, object$zeta)
}

nu.zicmp <- function(object, ...)
{
	exp(object$S %*% object$gamma)
}

sdev.zicmp <- function(object, ...)
{
	sqrt(diag(object$V))
}

# TBD: Update
chisq.zicmp <- function(object, ...)
{
	LRT(object$predictors, object$response, object$glm_coefficients, object$coef, object$nu, object$max)$teststat[1,1]
}

# TBD: Update
pval.zicmp <- function(object, ...)
{
	LRT(object$predictors, object$response, object$glm_coefficients, object$coef, object$nu, object$max)$pvalue
}

# TBD: Update
leverage.zicmp <- function(object, ...)
{
	CMPLeverage(object$predictors, object$response, object$coef, object$nu, object$max)
}

# TBD: Update
deviance.zicmp <- function(object, ...)
{
	CMPDeviance(object$predictors, object$response, object$coef, object$nu, leverage.cmp(object), object$max)
}

residuals.zicmp <- function(object, type = c("raw", "quantile"), ...)
{
	lambda.hat <- exp(object$X %*% object$beta)
	nu.hat <- exp(object$S %*% object$gamma)
	p.hat <- plogis(object$W %*% object$zeta)
	y.hat <- expected.y(lambda.hat, nu.hat, p.hat, object$max)

	type <- match.arg(type)
	if (type == "raw") {
		res <- object$y - y.hat
	} else if (type == "quantile") {
		res <- rqres.zicmp(object$y, lambda.hat, nu.hat, p.hat)
	} else {
		stop("Unsupported residual type")
	}

	return(as.numeric(res))
}

# TBD: Update
predict.zicmp <- function(object, ...)
{
	newdata = list(...)[["newdata"]]
	return(constantCMPfitsandresids(object$coef, object$nu, newdata[,object$x_names])$fit)
}

# TBD: Update
parametric_bootstrap.zicmp <- function(object, ...)
{
	n = list(...)[["n"]]
	if (is.null(n)) {
		n = 1000
	}
	bootstrap_results = as.data.frame(CMPParamBoot(x=object$predictors, object$glm_coefficients, betahat=object$coef, nuhat=object$nu, n=n)$CMPresult)
	names(bootstrap_results) = c("(Intercept)",object$x_names,"nu",recursive=TRUE)
	return(bootstrap_results)
}
