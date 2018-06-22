fit.zicmp.reg <- function(y, X, S, W, beta.init, gamma.init, zeta.init)
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
		logz <- z_hybrid(lambda, nu, take_log = TRUE)
		t(u) %*% log(p*exp(logz) + (1-p)) + t(1 - u) %*% (log(1-p) +
			y*log(lambda) - nu*lgamma(y+1)) - sum(logz)
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

	# TBD: Clean this up
	if (FALSE) {
		FIM <- fim.zicmp.reg(X, S, W, theta.hat$beta, theta.hat$gamma,
			theta.hat$zeta)
	} else {
		H <- optimHess(res$par, loglik, control = optim.control)
		FIM <- -H
	}
	colnames(FIM) <- rownames(FIM) <- c(
		sprintf("X:%s", colnames(X)),
		sprintf("S:%s", colnames(S)),
		sprintf("W:%s", colnames(W)))
	V <- solve(FIM)

	lambda.hat <- exp(X %*% theta.hat$beta)
	nu.hat <- exp(S %*% theta.hat$gamma)
	p.hat <- plogis(W %*% theta.hat$zeta)
	mu.hat <- expected.y(lambda.hat, nu.hat, p.hat)
	mse <- mean( (y - mu.hat)^2 )

	loglik <- res$value
	elapsed.sec <- as.numeric(Sys.time() - start, type = "sec")

	res <- list(theta.hat = theta.hat, V = V, H = H, FIM = FIM, opt.res = res,
		elapsed.sec = elapsed.sec, loglik = loglik, n = n)
	return(res)
}

fit.cmp.reg <- function(y, X, S, beta.init, gamma.init)
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
		logz <- z_hybrid(lambda, nu, take_log = TRUE)
		sum(y*log(lambda) - nu*lgamma(y+1) - logz)
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

	# TBD: Clean this up
	if (FALSE) {
		W <- matrix(1, n, 1)
		FIM.full <- fim.zicmp.reg(X, S, W = W, theta.hat$beta, theta.hat$gamma,
			zeta = -Inf)
		FIM <- FIM.full[1:(d1+d2), 1:(d1+d2)]
	} else {
		H <- optimHess(res$par, loglik, control = optim.control)
		FIM <- -H
	}
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
	mu.hat <- expected.y(lambda.hat, nu.hat, p.hat)
	mse <- mean( (y - mu.hat)^2 )

	loglik <- res$value
	elapsed.sec <- as.numeric(Sys.time() - start, type = "sec")

	res <- list(theta.hat = theta.hat, V = V, FIM = FIM, opt.res = res,
		elapsed.sec = elapsed.sec, loglik = loglik, n = n)
	return(res)
}

fit.zip.reg <- function(y, X, W, beta.init, zeta.init)
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

	# TBD: Clean this up
	if (FALSE) {
		FIM <- fim.zicmp.reg(X, S = matrix(1, n, 1), W, theta.hat$beta,
			gamma = 0, theta.hat$zeta)
	} else {
		H <- optimHess(res$par, loglik, control = optim.control)
		FIM <- -H
	}

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
	mu.hat <- expected.y(lambda.hat, nu.hat, p.hat)
	mse <- mean( (y - mu.hat)^2 )

	loglik <- res$value
	elapsed.sec <- as.numeric(Sys.time() - start, type = "sec")

	res <- list(theta.hat = theta.hat, V = V, FIM = FIM, opt.res = res,
		elapsed.sec = elapsed.sec, loglik = loglik, n = n)
	return(res)
}

expected.y <- function(lambda, nu, p)
{
	(1-p) * cmp_expected_value(lambda, nu)
}

expected.y.reg <- function(X, S, W, beta, gamma, zeta)
{
	lambda <- exp(X %*% beta)
	nu <- exp(S %*% gamma)
	p <- plogis(W %*% zeta)
	expected.y(lambda, nu, p)
}
