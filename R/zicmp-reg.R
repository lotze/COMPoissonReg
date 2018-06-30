summary.zicmp <- function(object, ...)
{
	n <- nrow(object$X)
	d1 <- ncol(object$X)
	d2 <- ncol(object$S)
	d3 <- ncol(object$W)

	V <- vcov(object)
	est <- coef(object)
	se <- sdev(object)
	z.val <- est / se
	p.val <- 2*pnorm(-abs(z.val))
	qq <- length(est)

	DF <- data.frame(
		Estimate = round(est, 4),
		SE = round(se, 4),
		z.value = round(z.val, 4),
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
		se <- sqrt(t(J) %*% V %*% J)

		DF.lambda <- data.frame(
			Estimate = round(est, 4),
			SE = round(se, 4)
		)
		rownames(DF.lambda) <- "lambda"
	}

	if (is.intercept.only(object$S)) {
		est <- exp(object$gamma)
		J <- c(rep(0, d1), exp(object$gamma), rep(0, d3))
		se <- sqrt(t(J) %*% V %*% J)

		DF.nu <- data.frame(
			Estimate = round(est, 4),
			SE = round(se, 4)		)
		rownames(DF.nu) <- "nu"
	}

	if (is.intercept.only(object$W)) {
		est <- plogis(object$zeta)
		J <- c(rep(0, d1), rep(0, d2), dlogis(object$zeta))
		se <- sqrt(t(J) %*% V %*% J)

		DF.p <- data.frame(
			Estimate = round(est, 4),
			SE = round(se, 4)
		)
		rownames(DF.p) <- "p"
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

print.zicmp <- function(x, ...)
{
	printf("ZICMP coefficients\n")
	s <- summary(x)
	tt <- equitest(x)
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

sdev.zicmp <- function(object, use.fim = FALSE, ...)
{
	sqrt(diag(vcov(object, use.fim)))
}

vcov.zicmp <- function(object, use.fim = FALSE, ...)
{
	if (use.fim) {
		# Compute the covariance via Fisher Information Matrix
		n <- nrow(object$X)
		X <- object$X
		S <- object$S
		W <- object$W
		FIM <- fim.zicmp.reg(X, S, W, object$beta, object$gamma,
			object$zeta)
		V <- solve(FIM)
		rownames(V) <- colnames(V) <- names(coef(object))
	} else {
		# Compute the covariance via Hessian from optimizer
		V <- -solve(object$H)
	}

	return(V)
}

equitest.zicmp <- function(object, ...)
{
	if ("equitest" %in% names(object)) {
		return(object$equitest)
	}

	y <- object$y
	X <- object$X
	S <- object$S
	W <- object$W
	beta.hat <- object$beta
	gamma.hat <- object$gamma
	zeta.hat <- object$zeta
	n <- length(y)

	fit0.out <- fit.zip.reg(y, X, W, beta.init = beta.hat,
		zeta.init = zeta.hat)
	beta0.hat <- fit0.out$theta.hat$beta
	zeta0.hat <- fit0.out$theta.hat$zeta

	ff <- numeric(n)
	ff0 <- numeric(n)

	lambda.hat <- exp(X %*% beta.hat)
	nu.hat <- exp(S %*% gamma.hat)
	p.hat <- plogis(W %*% zeta.hat)

	lambda0.hat <- exp(X %*% beta0.hat)
	nu0.hat <- rep(1, n)
	p0.hat <- plogis(W %*% zeta0.hat)

	ff <- dzicmp(y, lambda.hat, nu.hat, p.hat)
	ff0 <- dzicmp(y, lambda0.hat, nu0.hat, p0.hat)

	X2 <- 2*(sum(log(ff)) - sum(log(ff0)))
	pvalue <- pchisq(X2, df = length(gamma.hat), lower.tail = FALSE)
	list(teststat = X2, pvalue = pvalue)
}

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
			dzicmp(y[i], lambda, nu, p, log = TRUE)
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
		res <- rqres.zicmp(object$y, lambda.hat, nu.hat, p.hat)
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
	y.hat <- expected.y(lambda.hat, nu.hat, p.hat)
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
			logger("Starting bootstrap rep %d\n", r)
		}

		# Generate bootstrap samples of the full dataset using MLE
		y.boot <- rzicmp(n, lambda.hat, nu.hat, p.hat)

		# Take each of the bootstrap samples, along with the x matrix, and fit model
		# to generate bootstrap estimates
		tryCatch({
			fit.boot <- fit.zicmp.reg(y.boot, object$X, object$S, object$W,
				object$beta.init, object$gamma.init, object$zeta.init)
			theta.boot[r,] <- unlist(fit.boot$theta.hat)
		},
		error = function(e) {
			# Do nothing now; emit a warning later
		})
	}

	cnt <- sum(rowSums(is.na(theta.boot)) > 0)
	if (cnt > 0) {
		warning(sprintf("%d out of %d bootstrap iterations failed", cnt, reps))
	}

	colnames(theta.boot) <- c(
		sprintf("X:%s", colnames(object$X)),
		sprintf("S:%s", colnames(object$S)),
		sprintf("W:%s", colnames(object$W))
	)

	return(theta.boot)
}
