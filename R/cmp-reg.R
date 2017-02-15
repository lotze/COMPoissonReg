glm.cmp.old <- function(formula, initial.est=NULL, nuinit=1, max=100, ...){
	initial_glm = glm(formula, family='poisson', ...)
	if (is.null(initial.est)) {
		initial.est = coef(initial_glm)
	}

	# Parse formula
	mf <- model.frame(formula, ...)
	y <- model.response(mf)
	X <- model.matrix(formula, mf)

	object_result = list()
	object_result$formula = formula
	object_result$response = y
	object_result$predictors = X
	object_result$max = max

	internal_result = ComputeBetasAndNuHat(object_result$predictors, object_result$response, betainit=initial.est, nuinit=nuinit, max=object_result$max)
	if(internal_result$convergence == 1) {stop(sprintf("Constant CMP estimates could not be determined.  Optimization scheme did not converge: ", internal_result))}

	num_pars = length(internal_result$par)
	object_result$glm_coefficients = coef(initial_glm)
	object_result$coefficients = internal_result$par[1:(num_pars -1)]
	object_result$nu = internal_result$par[num_pars]
	object_result$loglik <- -internal_result$objective
	object_result$convergence <- internal_result$convergence
	object_result$message <- internal_result$message

	attr(object_result, "class") <- c("cmp", attr(object_result, "class"))
	return(object_result)
}

ComputeBetasAndNuHat <- function(x, y, betainit, nuinit, max)
{
	# Uses nlminb to solve for the MLE estimates for betas and nu
	# Create -logL so that we take the minimum of this function (which equals the max of logL)
	minusloglike <- function(par) {
		beta <- par[1:length(betainit)]
		nu <- par[length(betainit)+1]
		z <- computez(exp(x %*% beta), nu, max)
		-sum((y * (x %*% beta)) - (nu * lgamma(y+1)) - log(z))
	}

	# Determine the MLEs
	BetaNuEst <- nlminb(start=c(betainit,nuinit),minusloglike,lower = c(rep(-Inf,length(betainit)),0), upper = c(rep(Inf,length(betainit)),Inf))

	if (BetaNuEst$convergence != 0) {
		# Determine the MLEs
		BetaNuEst <- optim(par=c(betainit,nuinit),fn=minusloglike,control=list(maxit=1000))#$par
	}

	return(BetaNuEst)
}

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
		z.value = round(z.val, 6),
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
		z.val <- est / se
		p.val <- 2*pnorm(-abs(z.val))

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
		J <- c(rep(0, d1), exp(object$gamma))
		se <- sqrt(t(J) %*% object$V %*% J)
		z.val <- est / se
		p.val <- 2*pnorm(-abs(z.val))

		DF.nu <- data.frame(
			Estimate = round(est, 4),
			SE = round(se, 4),
			z.value = round(z.val, 6),
			p.value = sprintf("%0.4g", p.val)
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
	cat("Fit for CMP model\n")
	s <- summary.cmp(x)
	print(s$DF)

	if (!is.null(s$DF.lambda) || !is.null(s$DF.nu)) {
		cat("--\n")
		cat("Estimates for non-regression parameters\n")
		print(rbind(s$DF.lambda, s$DF.nu))
	}

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
	y <- object$y
	X <- object$X
	S <- object$S
	lambda0 <- exp(X %*% object$beta.glm)
	lambda <- exp(X %*% object$beta)
	nu <- exp(S %*% object$gamma)

	z <- computez(lambda, nu, max)
	teststat <- -2 * sum(y*log(lambda0) - lgamma(y+1) - lambda0 -
		y*log(lambda) + nu*lgamma(y+1) + log(z))
	pvalue <- pchisq(teststat, df=1, lower.tail=FALSE)
	list(teststat = teststat, pvalue = pvalue)
}

leverage.cmp <- function(object, ...)
{
	x <- object$X
	betahat <- object$beta
	nuhat <- exp(object$S %*% object$gamma)
	max <- object$max

	# 1) to code the W matrix  (diagonal matrix with Var(Y_i) )
	W <- diag(weights(x,betahat,nuhat,max))

	#    and X matrix (in Appendix)
	E.y <- computez.prodj(exp(x %*% betahat),nuhat,max)/computez(exp(x %*% betahat),nuhat,max)
	E.logfacty <- computez.prodlogj(exp(x %*% betahat),nuhat,max)/computez(exp(x %*% betahat),nuhat,max)
	extravec <- (-log(factorial(y)) + E.logfacty)/(y - E.y)
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
		minuslogf <- function(par){
			beta <- par[1:length(betainit)]
			ll <- y[i]*(x[i,] %*% beta) - nuhat[i]*lgamma(y[i]+1) -
				log(computez(exp(x[i,] %*% beta), nuhat[i], max))
			return(-ll)
		}
		
		# Determine the MLEs
		BetaEstResult <- nlm(p=betainit, f=minuslogf)
		OptimalLogLi[i] <- -BetaEstResult$min
	}

	#### Compute exact deviances
	lambdahat <- exp(x %*% betahat)
	OptimalLogL.mu <- (y*log(lambdahat)) - (nuhat * log(factorial(y))) - log(computez(lambdahat,nuhat,max))
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
	} else {
		X <- object$X
	}

	nu <- exp(object$S %*% object$gamma)
	out <- constantCMPfitsandresids(object$beta, nu, X)
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
	boot.out <- matrix(0, nrow=reps, ncol=length(betahat)+length(object$gamma))

	# Take each of the 1000 sample results along with the x matrix, and run CMP
	# regression on it to generate new betas and nu
	for (i in 1:reps){
		if (i %% report.period == 0) {
			logger("Starting bootstrap rep %d\n", i)
		}
		Ystar[,i] <- rcmp(nrow(X), lambdahat, nuhat, max = max)
		res <- fit.cmp.reg(y = Ystar[,i], X = X, S = S,
			beta.init = poissonest, gamma.init = object$gamma, max = object$max)
		boot.out[i,] <- unlist(res$theta.hat)
	}

	colnames(boot.out) <- c(colnames(X), colnames(S), recursive=TRUE)
	return(boot.out)
}

constantCMPfitsandresids <- function(betahat, nuhat, x, y=0)
{
	# Determine estimated lambdahat, fit, and residuals
	lambdahat <- exp(x %*% betahat)
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
	CMPresids <- constantCMPfitsandresids(CMPbetas,CMPnu,x,y)$resid
	MSE <- mean(CMPresids^2)
	return(MSE)
}

