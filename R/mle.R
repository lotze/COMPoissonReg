fit.zicmp.reg <- function(y, X, S, W, beta.init, gamma.init, zeta.init, max)
{
	start <- Sys.time()
	u <- as.integer(y == 0)
	d1 <- ncol(X)
	d2 <- ncol(S)
	d3 <- ncol(W)
	n <- length(y)
	qq <- d1 + d2 + d3

	if (is.null(beta.init)) { beta.init <- rep(0, d1) }
	if (is.null(gamma.init)) { gamma.init <- rep(0, d2) }
	if (is.null(zeta.init)) { zeta.init <- rep(0, d3) }
	stopifnot(d1 == length(beta.init))
	stopifnot(d2 == length(gamma.init))
	stopifnot(d3 == length(zeta.init))

	tx <- function(par) {
		list(
			beta = par[1:d1],
			gamma = par[1:d2 + d1],
			zeta = par[1:d3 + d1 + d2]
		)
	}

	loglik <- function(par) {
		theta <- tx(par)
		lambda <- exp(X %*% theta$beta)
		nu <- exp(S %*% theta$gamma)
		p <- plogis(W %*% theta$zeta)
		logz <- computez(lambda, nu, max, log = TRUE, autoscale = TRUE)
		t(u) %*% log(p*exp(logz) + (1-p)) + t(1 - u) %*% (log(1-p) +
			y*log(lambda) - nu*lgamma(y+1)) - sum(logz)

		# browser()
		# logz_old <- computez(lambda, nu, max, log = TRUE)
		# t(u) %*% log(p*exp(logz_old) + (1-p)) + t(1 - u) %*% (log(1-p) +
		# 	y*log(lambda) - nu*lgamma(y+1)) - sum(logz_old)
	}

	optim.method <- getOption("COMPoissonReg.optim.method")
	optim.control <- getOption("COMPoissonReg.optim.control")

	optim.control$fnscale = -1
	par.init <- c(beta.init, gamma.init, zeta.init)
	res <- optim(par.init, loglik, method = optim.method,
		control = optim.control, hessian = TRUE)
	H <- res$hessian

	theta.hat <- list(
		beta = res$par[1:d1],
		gamma = res$par[1:d2 + d1],
		zeta = res$par[1:d3 + d1 + d2]
	)
	names(theta.hat$beta) <- sprintf("X:%s", colnames(X))
	names(theta.hat$gamma) <- sprintf("S:%s", colnames(S))
	names(theta.hat$zeta) <- sprintf("W:%s", colnames(W))

	FIM <- fim.zicmp.reg(X, S, W, theta.hat$beta, theta.hat$gamma,
		theta.hat$zeta, max = max)
	colnames(FIM) <- rownames(FIM) <- c(
		sprintf("X:%s", colnames(X)),
		sprintf("S:%s", colnames(S)),
		sprintf("W:%s", colnames(W)))
	V <- solve(FIM)

	lambda.hat <- exp(X %*% theta.hat$beta)
	nu.hat <- exp(S %*% theta.hat$gamma)
	p.hat <- plogis(W %*% theta.hat$zeta)
	mu.hat <- expected.y(lambda.hat, nu.hat, p.hat, max)
	mse <- mean( (y - mu.hat)^2 )

	loglik <- res$value
	elapsed.sec <- as.numeric(Sys.time() - start, type = "sec")

	res <- list(theta.hat = theta.hat, V = V, H = H, FIM = FIM, opt.res = res,
		elapsed.sec = elapsed.sec, loglik = loglik, n = n)
	return(res)
}

fit.cmp.reg <- function(y, X, S, beta.init, gamma.init, max)
{
	start <- Sys.time()
	u <- as.integer(y == 0)
	d1 <- ncol(X)
	d2 <- ncol(S)
	n <- length(y)
	qq <- d1 + d2

	if (is.null(beta.init)) { beta.init <- rep(0, d1) }
	if (is.null(gamma.init)) { gamma.init <- rep(0, d2) }
	stopifnot(d1 == length(beta.init))
	stopifnot(d2 == length(gamma.init))

	tx <- function(par) {
		list(
			beta = par[1:d1],
			gamma = par[1:d2 + d1]
		)
	}

	loglik <- function(par) {
		theta <- tx(par)
		lambda <- exp(X %*% theta$beta)
		nu <- exp(S %*% theta$gamma)
		
		if (TRUE) {
			logz <- computez(lambda, nu, max, log = TRUE, autoscale = FALSE)
			ans <- sum(y*log(lambda) - nu*lgamma(y+1) - logz)
		} else {
			# logscale <- max*log(lambda)
			# logscale <- rep(40,n)
			logz <- computez(lambda, nu, max, log = TRUE, autoscale = TRUE)
			ans <- sum(y*log(lambda) - nu*lgamma(y+1) - logz)
			# if (is.na(ans) || is.infinite(ans)) { browser() }
		}
print(ans)
		return(ans)

		# browser()		
		# logz_old <- computez(lambda, nu, max, log = TRUE)
		# sum(y*log(lambda) - nu*lgamma(y+1) - logz_old$ans)
	}

	optim.method <- getOption("COMPoissonReg.optim.method")
	optim.control <- getOption("COMPoissonReg.optim.control")
	optim.control$fnscale = -1
	par.init <- c(beta.init, gamma.init)

	res <- optim(par.init, loglik, method = optim.method,
		control = optim.control)

	theta.hat <- list(
		beta = res$par[1:d1],
		gamma = res$par[1:d2 + d1]
	)
	names(theta.hat$beta) <- sprintf("X:%s", colnames(X))
	names(theta.hat$gamma) <- sprintf("S:%s", colnames(S))

	W <- matrix(1, n, 1)
	FIM.full <- fim.zicmp.reg(X, S, W = W, theta.hat$beta, theta.hat$gamma,
		zeta = -Inf, max = max)
	FIM <- FIM.full[1:(d1+d2), 1:(d1+d2)]
	colnames(FIM) <- rownames(FIM) <- c(
		sprintf("X:%s", colnames(X)),
		sprintf("S:%s", colnames(S)))

	V <- tryCatch({
		solve(FIM)
	}, error = function(e) {
		warning("FIM could not be inverted. Trying Hessian to estimate variance")
		H <- optimHess(res$par, loglik, control = optim.control)
		solve(-H)
	})

	lambda.hat <- exp(X %*% theta.hat$beta)
	nu.hat <- exp(S %*% theta.hat$gamma)
	p.hat <- 0
	mu.hat <- expected.y(lambda.hat, nu.hat, p.hat, max)
	mse <- mean( (y - mu.hat)^2 )

	loglik <- res$value
	elapsed.sec <- as.numeric(Sys.time() - start, type = "sec")

	res <- list(theta.hat = theta.hat, V = V, FIM = FIM, opt.res = res,
		elapsed.sec = elapsed.sec, loglik = loglik, n = n)
	return(res)
}

fit.zip.reg <- function(y, X, W, beta.init, zeta.init, max)
{
	start <- Sys.time()
	u <- as.integer(y == 0)
	d1 <- ncol(X)
	d3 <- ncol(W)
	n <- length(y)
	qq <- d1 + d3

	if (is.null(beta.init)) { beta.init <- rep(0, d1) }
	if (is.null(zeta.init)) { zeta.init <- rep(0, d3) }
	stopifnot(d1 == length(beta.init))
	stopifnot(d3 == length(zeta.init))

	tx <- function(par) {
		list(
			beta = par[1:d1],
			zeta = par[1:d3 + d1]
		)
	}

	loglik <- function(par) {
		theta <- tx(par)
		lambda <- exp(X %*% theta$beta)
		p <- plogis(W %*% theta$zeta)
		sum(u*log(p + (1-p)*exp(-lambda)) + (1-u)*log(1-p) +
			(1-u)*(y*log(lambda) - lambda - lgamma(y+1)))
	}

	optim.method <- getOption("COMPoissonReg.optim.method")
	optim.control <- getOption("COMPoissonReg.optim.control")

	optim.control$fnscale = -1
	par.init <- c(beta.init, zeta.init)
	res <- optim(par.init, loglik, method = optim.method,
		control = optim.control)

	theta.hat <- list(
		beta = res$par[1:d1],
		zeta = res$par[1:d3 + d1]
	)

	FIM <- fim.zicmp.reg(X, S = matrix(1, n, 1), W, theta.hat$beta,
		gamma = 0, theta.hat$zeta, max = max)

	V <- tryCatch({
		solve(FIM)
	}, error = function(e) {
		warning("FIM could not be inverted. Trying Hessian to estimate variance")
		H <- optimHess(res$par, loglik, control = optim.control)
		solve(-H)
	})

	lambda.hat <- exp(X %*% theta.hat$beta)
	nu.hat <- rep(0, n)
	p.hat <- plogis(W %*% theta.hat$zeta)
	mu.hat <- expected.y(lambda.hat, nu.hat, p.hat, max)
	mse <- mean( (y - mu.hat)^2 )

	loglik <- res$value
	elapsed.sec <- as.numeric(Sys.time() - start, type = "sec")

	res <- list(theta.hat = theta.hat, V = V, FIM = FIM, opt.res = res,
		elapsed.sec = elapsed.sec, loglik = loglik, n = n)
	return(res)
}

expected.y <- function(lambda, nu, p, max)
{
	z <- computez(lambda, nu, max=max)
	dzdlambda <- computez.prodj(lambda, nu, max=max)/lambda
	dlogzdlambda <- dzdlambda / z
	(1-p) * lambda * dlogzdlambda	
}

expected.y.reg <- function(X, S, W, beta, gamma, zeta, max)
{
	lambda <- exp(X %*% beta)
	nu <- exp(S %*% gamma)
	p <- plogis(W %*% zeta)
	expected.y(lambda, nu, p, max)
}
