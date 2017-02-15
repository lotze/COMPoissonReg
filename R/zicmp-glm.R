glm.zicmp <- function(formula.lambda, formula.nu = NULL, formula.p = NULL,
	beta.init = NULL, gamma.init = NULL, zeta.init = NULL, max = 100, ...)
{
	# Parse formula.lambda. This one should have the response.
	mf <- model.frame(formula.lambda, ...)
	y <- model.response(mf)
	X <- model.matrix(formula.lambda, mf)
	d1 <- ncol(X)

	if (is.null(formula.nu)) { formula.nu <- ~ 1 }
	if (is.null(formula.p)) { formula.p <- ~ 1 }
	
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

	res$beta.init <- beta.init
	res$gamma.init <- gamma.init
	res$zeta.init <- zeta.init

	fit.out <- fit.zicmp.reg(res$y, res$X, res$S, res$W, beta.init = beta.init,
		gamma.init = gamma.init, zeta.init = zeta.init, max = res$max)

	res$beta.glm <- coef(initial.glm)
	res$beta <- fit.out$theta.hat$beta
	res$gamma <- fit.out$theta.hat$gamma
	res$zeta <- fit.out$theta.hat$zeta
	res$FIM <- fit.out$FIM
	res$V <- fit.out$V
	res$loglik <- fit.out$loglik
	res$opt.res <- fit.out$opt.res
	res$elapsed.sec <- fit.out$elapsed.sec

	attr(res, "class") <- c("zicmp", attr(res, "class"))
	return(res)
}

summary.zicmp <- function(object, ...)
{
	n <- nrow(object$X)
	d1 <- ncol(object$X)
	d2 <- ncol(object$S)
	d3 <- ncol(object$W)

	est <- coef(object)
	se <- sdev(object)
	z.val <- est / se
	p.val <- 2*pnorm(-abs(z.val))
	qq <- length(est)

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

	# If X, S, or W are intercept only, compute results for non-regression parameters
	DF.lambda <- NULL
	DF.nu <- NULL
	DF.p <- NULL
	X.int <- matrix(1, n, 1)

	is.intercept.only <- function(A, eps = 1e-12) {
		prod(dim(A) == c(n,1)) & norm(A - 1, type = "F") < eps
	}

	if (is.intercept.only(object$X)) {
		est <- exp(object$beta)
		J <- c(exp(object$beta), rep(0, d2), rep(0, d3))
		se <- sqrt(t(J) %*% object$V %*% J)
		z.val <- est / se
		p.val <- 2*(1 - pnorm(abs(z.val)))

		DF.lambda <- data.frame(
			Estimate = round(est, 4),
			SE = round(se, 4),
			z.value = round(z.val, 6),
			p.value = sprintf("%0.4g", p.val)
		)
		rownames(DF.lambda) <- "lambda"
	}

	if (is.intercept.only(object$S)) {
		est <- exp(object$gamma)
		J <- c(rep(0, d1), exp(object$gamma), rep(0, d3))
		se <- sqrt(t(J) %*% object$V %*% J)
		z.val <- est / se
		p.val <- 2*(1 - pnorm(abs(z.val)))

		DF.nu <- data.frame(
			Estimate = round(est, 4),
			SE = round(se, 4),
			z.value = round(z.val, 6),
			p.value = sprintf("%0.4g", p.val)
		)
		rownames(DF.nu) <- "nu"
	}

	if (is.intercept.only(object$W)) {
		est <- plogis(object$zeta)
		J <- c(rep(0, d1), rep(0, d2), dlogis(object$zeta))
		se <- sqrt(t(J) %*% object$V %*% J)
		z.val <- est / se
		p.val <- 2*(1 - pnorm(abs(z.val)))

		DF.p <- data.frame(
			Estimate = round(est, 4),
			SE = round(se, 4),
			z.value = round(z.val, 6),
			p.value = sprintf("%0.4g", p.val)
		)
		rownames(DF.p) <- "p"
	}

	list(DF = DF, DF.lambda = DF.lambda, DF.nu = DF.nu, DF.p = DF.p,
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

	if (!is.null(s$DF.lambda) || !is.null(s$DF.nu) || !is.null(s$DF.p)) {
		cat("--\n")
		cat("Estimates for non-regression parameters\n")
		print(rbind(s$DF.lambda, s$DF.nu, s$DF.p))
	}

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

logLik.zicmp <- function(object, ...)
{
	object$loglik
}

AIC.zicmp <- function(object, ..., k = 2)
{
	-2*object$loglik + 2*length(coef(object))
}

BIC.zicmp <- function(object, ...)
{
	n <- length(object$y)
	-2*object$loglik + log(n)*length(coef(object))
}

coef.zicmp <- function(object, ...)
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

equitest.zicmp <- function(object, ...)
{
	fit0.out <- fit.zip.reg(object$y, object$X, object$W, beta.init = object$beta,
		zeta.init = object$zeta, max = object$max)
	res <- LRT.zicmp(object$y, object$X, object$S, object$W, object$beta, object$gamma,
		object$zeta, fit0.out$theta.hat$beta, fit0.out$theta.hat$zeta, object$max)
	list(teststat = res$stat, pvalue = res$pvalue)
}

# TBD: discuss this. How to define leverage for ZICMP?
leverage.zicmp <- function(object, ...)
{
	stop("This function is not yet implemented")
}

deviance.zicmp <- function(object, ...)
{
	y <- object$y
	X <- object$X
	S <- object$S
	W <- object$W
	n <- length(y)
	d1 <- ncol(X)
	d2 <- ncol(S)
	d3 <- ncol(W)

	par.hat <- c(object$beta, object$gamma, object$zeta)
	par.init <- par.hat
	ll <- numeric(n)
	ll.star <- numeric(n)

	for(i in 1:n) {
		loglik <- function(par){
			lambda <- exp(X[i,] %*% par[1:d1])
			nu <- exp(S[i,] %*% par[1:d2 + d1])
			p <- plogis(W[i,] %*% par[1:d3 + d1 + d2])
			dzicmp(y[i], lambda, nu, p, max = object$max, log = TRUE)
		}

		# Maximize loglik for ith obs
		res <- optim(par.init, loglik, control = list(fnscale = -1))
		ll.star[i] <- res$value

		# loglik maximized over all the obs, evaluated for ith obs
		ll[i] <- loglik(par.hat)
	}

	dev <- -2*(ll - ll.star)
	leverage <- leverage(object)
	cmpdev <- dev / sqrt(1 - leverage)
	return(cmpdev)
}

residuals.zicmp <- function(object, type = c("raw", "quantile"), ...)
{
	lambda.hat <- exp(object$X %*% object$beta)
	nu.hat <- exp(object$S %*% object$gamma)
	p.hat <- plogis(object$W %*% object$zeta)
	y.hat <- predict.zicmp(object)

	type <- match.arg(type)
	if (type == "raw") {
		res <- object$y - y.hat
	} else if (type == "quantile") {
		res <- rqres.zicmp(object$y, lambda.hat, nu.hat, p.hat, object$max)
	} else {
		stop("Unsupported residual type")
	}

	return(as.numeric(res))
}

predict.zicmp <- function(object, newdata = NULL, ...)
{
	if (!is.null(newdata)) {
		# If any of the original models had an intercept added via model.matrix, they
		# will have an "(Intercept)" column. Let's add an "(Intercept)" to newdata
		# in case the user didn't make one.
		newdata$'(Intercept)' <- 1

		X <- as.matrix(newdata[,colnames(object$X)])
		S <- as.matrix(newdata[,colnames(object$S)])
		W <- as.matrix(newdata[,colnames(object$W)])
	} else {
		X <- object$X
		S <- object$S
		W <- object$W
	}

	lambda.hat <- exp(X %*% object$beta)
	nu.hat <- exp(S %*% object$gamma)
	p.hat <- plogis(W %*% object$zeta)
	y.hat <- expected.y(lambda.hat, nu.hat, p.hat, object$max)
	return(y.hat)
}

parametric_bootstrap.zicmp <- function(object, reps = 1000, report.period = reps+1, ...)
{
	n <- length(object$y)
	qq <- length(object$beta) + length(object$gamma) + length(object$zeta)
	theta.boot <- matrix(NA, reps, qq)

	lambda.hat <- exp(object$X %*% object$beta)
	nu.hat <- exp(object$S %*% object$gamma)
	p.hat <- plogis(object$W %*% object$zeta)

	for (r in 1:reps) {
		if (r %% report.period == 0) {
			logger("Starting boostrap rep %d\n", r)
		}

		# Generate bootstrap samples of the full dataset using MLE
		y.boot <- rzicmp(n, lambda.hat, nu.hat, p.hat, max = object$max)

		# Take each of the bootstrap samples, along with the x matrix, and fit model
		# to generate bootstrap estimates
		tryCatch({
			fit.boot <- fit.zicmp.reg(y.boot, object$X, object$S, object$W,
				object$beta.init, object$gamma.init, object$zeta.init, object$max)
			theta.boot[r,] <- unlist(fit.boot$theta.hat)
		},
		error = function(e) {
			print(e)
			theta.boot[r,] <- NA
		})
	}

	colnames(theta.boot) <- c(
		sprintf("X:%s", colnames(object$X)),
		sprintf("S:%s", colnames(object$S)),
		sprintf("W:%s", colnames(object$W))
	)

	return(theta.boot)
}

