summary.cmp <- function(object, ...)
{
	n <- nrow(object$X)
	d1 <- ncol(object$X)
	d2 <- ncol(object$S)

	est <- coef(object)
	se <- sdev(object)
	z.val <- est / se
	p.val <- 2*pnorm(-abs(z.val))

	DF <- data.frame(
		Estimate = round(est, 4),
		SE = round(se, 4),
		z.value = round(z.val, 4),
		p.value = sprintf("%0.4g", p.val)
	)
	rownames(DF) <- c(sprintf("X:%s", colnames(object$X)),
		sprintf("S:%s", colnames(object$S)))

	# If X, S, or W are intercept only, compute results for non-regression parameters
	DF.lambda <- NULL
	DF.nu <- NULL
	X.int <- matrix(1, n, 1)

	is.intercept.only <- function(A, eps = 1e-12) {
		prod(dim(A) == c(n,1)) & norm(A - 1, type = "F") < eps
	}

	if (is.intercept.only(object$X)) {
		est <- exp(object$beta)
		J <- c(exp(object$beta), rep(0, d2))
		se <- sqrt(t(J) %*% object$V %*% J)

		DF.lambda <- data.frame(
			Estimate = round(est, 4),
			SE = round(se, 4)
		)
		rownames(DF.lambda) <- "lambda"
	}

	if (is.intercept.only(object$S)) {
		est <- exp(object$gamma)
		J <- c(rep(0, d1), exp(object$gamma))
		se <- sqrt(t(J) %*% object$V %*% J)

		DF.nu <- data.frame(
			Estimate = round(est, 4),
			SE = round(se, 4)
		)
		rownames(DF.nu) <- "nu"
	}

	list(DF = DF, DF.lambda = DF.lambda, DF.nu = DF.nu,
		 n = n,
		 loglik = logLik(object),
		 aic = AIC(object),
		 bic = BIC(object),
		 opt.message = object$message,
		 opt.convergence = object$convergence,
		 elapsed.sec = object$elapsed.sec
	)
}

print.cmp <- function(x, ...)
{
	cat("Fit for CMP coefficients\n")
	s <- summary.cmp(x)
	tt <- equitest(x)
	print(s$DF)

	if (!is.null(s$DF.lambda) || !is.null(s$DF.nu)) {
		cat("--\n")
		cat("Transformed intercept-only parameters\n")
		print(rbind(s$DF.lambda, s$DF.nu))
	}
	cat("--\n")
	cat("Chi-squared test for equidispersion\n")
	cat(sprintf("X^2 = %0.4f, df = 1, ", tt$teststat))
	cat(sprintf("p-value = %0.4e\n", tt$pvalue))
	cat("--\n")
	cat(sprintf("Elapsed Sec: %0.2f   ", s$elapsed.sec))
	cat(sprintf("Sample size: %d\n", s$n))
	cat(sprintf("LogLik: %0.4f   ", s$loglik))
	cat(sprintf("AIC: %0.4f   ", s$aic))
	cat(sprintf("BIC: %0.4f   ", s$bic))
	cat("\n")
	cat(sprintf("Converged status: %d   ", s$opt.convergence))
	cat(sprintf("Message: %s\n", s$opt.message))
}

logLik.cmp <- function(object, ...)
{
	object$loglik
}

AIC.cmp <- function(object, ..., k=2)
{
	-2*object$loglik + 2*length(coef(object))
}

BIC.cmp <- function(object, ...)
{
	n <- length(object$y)
	-2*object$loglik + log(n)*length(coef(object))
}

coef.cmp <- function(object, ...)
{
	c(object$beta, object$gamma)
}

nu.cmp <- function(object, ...)
{
	exp(object$S %*% object$gamma)
}

sdev.cmp <- function(object, ...)
{
	sqrt(diag(object$V))
}

equitest.cmp <- function(object, ...)
{
	if ("equitest" %in% names(object)) {
		return(object$equitest)
	}

	y <- object$y
	X <- object$X
	S <- object$S
	lambda0 <- exp(X %*% object$beta.glm)
	lambda <- exp(X %*% object$beta)
	nu <- exp(S %*% object$gamma)
	max <- object$max

	z <- computez(lambda, nu, max)
	teststat <- -2 * sum(y*log(lambda0) - lgamma(y+1) - lambda0 -
		y*log(lambda) + nu*lgamma(y+1) + log(z))
	pvalue <- pchisq(teststat, df=1, lower.tail=FALSE)
	list(teststat = teststat, pvalue = pvalue)
}

leverage.cmp <- function(object, ...)
{
	y <- object$y
	x <- object$X
	betahat <- object$beta
	nuhat <- exp(object$S %*% object$gamma)
	max <- object$max

	# 1) to code the W matrix  (diagonal matrix with Var(Y_i) )
	W <- diag(weights(x,betahat,nuhat,max))

	#    and X matrix (in Appendix)
	E.y <- computez.prodj(exp(x %*% betahat),nuhat,max)/computez(exp(x %*% betahat),nuhat,max)
	E.logfacty <- computez.prodlogj(exp(x %*% betahat),nuhat,max)/computez(exp(x %*% betahat),nuhat,max)
	extravec <- (-lgamma(y+1) + E.logfacty)/(y - E.y)
	curlyX.mat <- cbind(x,extravec)
	
	# 2) to compute H using eq (12)  on p. 11
	H1 <- t(curlyX.mat) %*% sqrt(W)
	H2 <- solve(t(curlyX.mat) %*% W %*% curlyX.mat)
	H <- t(H1) %*% H2 %*% H1
	diagH <- diag(H)
	return(diagH)
}

deviance.cmp <- function(object, ...)
{
	# Compute the COM-Poisson deviances exactly
	y <- object$y
	x <- object$X
	betahat <- object$beta
	nuhat <- exp(object$S %*% object$gamma)
	max <- object$max
	leverage <- leverage.cmp(object)

	#### Compute optimal log likelihood value for given nu-hat value
	betainit <- betahat
	OptimalLogLi <- rep(0,length(y))
	iterct <- rep(0,length(y))

	for (i in 1:length(y)){
		# Create -logL = -logf (because considering single observation) so that we
		# take the minimum of this function (which equals the max of logL)
		logf <- function(par){
			beta <- par[1:length(betainit)]
			lambda <- exp(x[i,] %*% beta)
			ll <- y[i]*log(lambda) - nuhat[i]*lgamma(y[i]+1) -
				log(computez(lambda, nuhat[i], max))
			return(ll)
		}
		
		# Determine the MLEs
		# Using optim rather than nlm because optim handles -Inf more gracefully, if encountered
		BetaEstResult <- optim(betainit, logf, method = "L-BFGS-B", control = list(fnscale = -1))
		OptimalLogLi[i] <- BetaEstResult$value
	}

	#### Compute exact deviances
	lambdahat <- exp(x %*% betahat)
	OptimalLogL.mu <- (y*log(lambdahat)) - (nuhat * lgamma(y+1)) - log(computez(lambdahat,nuhat,max))
	OptimalLogL.y <- OptimalLogLi
	d <- -2*(OptimalLogL.mu - OptimalLogL.y)
	cmpdev <- d/(sqrt(1-leverage))
	return(cmpdev)
}

residuals.cmp <- function(object, type = c("raw", "quantile"), ...)
{
	X <- object$X
	S <- object$S
	lambda.hat <- exp(X %*% object$beta)
	nu.hat <- exp(S %*% object$gamma)
	y.hat <- predict.cmp(object, newdata = object$X)

	type <- match.arg(type)
	if (type == "raw") {
		res <- object$y - y.hat
	} else if (type == "quantile") {
		set.seed(1234)
		res <- rqres.cmp(object$y, lambda = lambda.hat, nu = nu.hat, max = object$max)
	} else {
		stop("Unsupported residual type")
	}

	return(as.numeric(res))
}

predict.cmp <- function(object, newdata = NULL, ...)
{
	if (!is.null(newdata)) {
		# If any of the original models had an intercept added via model.matrix, they
		# will have an "(Intercept)" column. Let's add an "(Intercept)" to newdata
		# in case the user didn't make one.
		newdata <- as.data.frame(newdata)
        newdata$'(Intercept)' <- 1
		X <- as.matrix(newdata[,colnames(object$X)])
		S <- as.matrix(newdata[,colnames(object$S)])
	} else {
		X <- object$X
		S <- object$S
	}

	lambda <- exp(X %*% object$beta)
	nu <- exp(S %*% object$gamma)
	out <- constantCMPfitsandresids(lambda, nu)
	y.hat <- out$fit
	return(y.hat)
}

parametric_bootstrap.cmp <- function(object, reps = 1000, report.period = reps+1, ...)
{
	n <- reps
	X <- object$X
	S <- object$S
	poissonest <- object$beta.glm
	betahat <- object$beta
	lambdahat <- exp(X %*% betahat)
	nuhat <- exp(object$S %*% object$gamma)
	max <- object$max

	# Generate 1000 samples, using betahat and nuhat from full dataset
	Ystar <- matrix(0, nrow = nrow(X), ncol = reps)
	boot.out <- matrix(NA, nrow=reps, ncol=length(betahat)+length(object$gamma))

	# Take each of the 1000 sample results along with the x matrix, and run CMP
	# regression on it to generate new betas and nu
	for (i in 1:reps){
		if (i %% report.period == 0) {
			logger("Starting bootstrap rep %d\n", i)
		}
		Ystar[,i] <- rcmp(nrow(X), lambdahat, nuhat, max = max)
		tryCatch({
			res <- fit.cmp.reg(y = Ystar[,i], X = X, S = S,
			   beta.init = poissonest, gamma.init = object$gamma, max = object$max)
			boot.out[i,] <- unlist(res$theta.hat)
		}, error = function(e) {
			# Do nothing now; emit a warning later
		})
	}

	cnt <- sum(rowSums(is.na(boot.out)) > 0)
	if (cnt > 0) {
		warning(sprintf("%d out of %d bootstrap iterations failed", cnt, reps))
	}

	colnames(boot.out) <- c(colnames(X), colnames(S), recursive=TRUE)
	return(boot.out)
}

constantCMPfitsandresids <- function(lambdahat, nuhat, y=0)
{
	# Determine estimated lambdahat, fit, and residuals
	lambdahat <- as.numeric(lambdahat)
	nuhat <- as.numeric(nuhat)
	fit <- lambdahat^(1/nuhat) - ((nuhat - 1)/(2*nuhat))
	resid <- y - fit
	list(fit = fit, resid = resid)
}

weights <- function(x,beta,nu,max)
{
	# Compute the parts that comprise the weight functions
	w1 <- computez.prodj2(exp(x %*% beta),nu,max)
	w2 <- computez.prodj(exp(x %*% beta),nu,max)
	w3 <- computez(exp(x %*% beta),nu,max)

	Ey2 <- w1/w3
	E2y <- (w2/w3)^2

	# Determine the weights
	weight <- Ey2 - E2y
	return(weight)
}

CMP.MSE <- function(CMPbetas,CMPnu,x,y)
{
	lambda <- exp(x %*% CMPbetas)
	CMPresids <- constantCMPfitsandresids(lambda,CMPnu,y)$resid
	MSE <- mean(CMPresids^2)
	return(MSE)
}

