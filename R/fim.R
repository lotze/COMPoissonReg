fim.cmp <- function(lambda, nu)
{
	f <- function(theta) {
		z_hybrid(theta[1], theta[2], take_log = TRUE)
	}
	H <- hessian(func = f, x = c(lambda, nu))
	H[1,1] <- H[1,1] + cmp_expected_value(lambda, nu) / lambda^2
	return(H)
}

fim.zicmp <- function(lambda, nu, p)
{
	stopifnot(length(lambda) == 1)
	stopifnot(length(nu) == 1)
	stopifnot(length(p) == 1)
	
	z <- z_hybrid(lambda, nu)
	meany <- zicmp_expected_value(lambda, nu, p)
	f <- function(theta) {
		z_hybrid(theta[1], theta[2], take_log = TRUE)
	}

	FIM <- matrix(NA, 3, 3)
	colnames(FIM) <- c("lambda", "nu", "p")
	rownames(FIM) <- c("lambda", "nu", "p")

	# The gradient and hessian of the z-like functions all have a closed
	# form, but require some work to compute reliably. We have some code to
	# compute the z-function itself somewhat reliably, so it's safer
	# to use numerical derivatives here. We use forward derivatives
	# to avoid issues at the boundary when lambda or nu is near 0.
	Dlogz <- grad.fwd(f, c(lambda, nu))
	Hlogz <- hess.fwd(f, c(lambda, nu))

	dlogzdlambda <- Dlogz[1]
	dlogzdnu <- Dlogz[2]
	d2logzdlambda2 <- Hlogz[1,1]
	d2logzdnu2 <- Hlogz[2,2]
	d2logzdlambdadnu <- Hlogz[1,2]

	# FIM[lambda, lambda]
	FIM[1,1] <- (1-p)*d2logzdlambda2 - p*(1-p)*dlogzdlambda^2 / (p*(z-1)+1) + meany / lambda^2

	# FIM[nu, nu]
	FIM[2,2] <- (1-p)*d2logzdnu2 - p*(1-p)*dlogzdnu^2 / (p*(z-1)+1)

	# FIM[p, p]
	FIM[3,3] <- (1/z) * (z-1)^2 / (p*(z-1) + 1) + (1 - 1/z) / (1-p)

	# FIM[lambda, nu]
	FIM[1,2] <- (1-p)*d2logzdlambdadnu - p*(1-p)*dlogzdnu*dlogzdlambda / (p*(z-1)+1)
	FIM[2,1] <- FIM[1,2]

	# FIM[lambda, p]
	FIM[1,3] <- -1/(p*(z-1) + 1) * dlogzdlambda
	FIM[3,1] <- FIM[1,3]

	# FIM[nu, p]
	FIM[2,3] <- -1/(p*(z-1) + 1) * dlogzdnu
	FIM[3,2] <- FIM[2,3]

	return(FIM)
}

fim.zicmp.reg <- function(X, S, W, beta, gamma, zeta)
{
	n <- nrow(X)
	qq <- ncol(X) + ncol(S) + ncol(W)
	lambda <- as.numeric(exp(X %*% beta))
	nu <- as.numeric(exp(S %*% gamma))
	p <- as.numeric(plogis(W %*% zeta))

	FIM <- matrix(0, qq, qq)
	colnames(FIM) <- c(
		sprintf("beta%d", 1:length(beta)),
		sprintf("gamma%d", 1:length(gamma)),
		sprintf("zeta%d", 1:length(zeta)))
	rownames(FIM) <- colnames(FIM)

	# Compute the FIM with respect to each (lambda, nu, p), then
	# transform to the FIM of (beta, gamma, zeta)
	FIM.one.list <- list()
	for (i in 1:n) {
		FIM.one.list[[i]] <- fim.zicmp(lambda[i], nu[i], p[i])
	}

	idx.beta <- 1:length(beta)
	idx.gamma <- 1:length(gamma) + length(beta)
	idx.zeta <- 1:length(zeta) + length(gamma) + length(beta)

	# FIM[beta, beta]
	D <- unlist(Map(function(x) { x[1,1] }, FIM.one.list))
	FIM[idx.beta, idx.beta] <- t(X) %*% ((D * lambda^2) * X)

	# FIM[gamma, gamma]
	D <- unlist(Map(function(x) { x[2,2] }, FIM.one.list))
	FIM[idx.gamma, idx.gamma] <- t(S) %*% ((D * nu^2) * S)

	# FIM[zeta, zeta]
	D <- unlist(Map(function(x) { x[3,3] }, FIM.one.list))
	FIM[idx.zeta, idx.zeta] <- t(W) %*% ((D * p^2*(1-p)^2) * W)

	# FIM[beta, gamma]
	D <- unlist(Map(function(x) { x[1,2] }, FIM.one.list))
	FIM[idx.beta, idx.gamma] <- t(X) %*% (D * lambda * nu * S)
	FIM[idx.gamma, idx.beta] <- t(FIM[idx.beta, idx.gamma])

	# FIM[zeta, gamma]
	D <- unlist(Map(function(x) { x[3,2] }, FIM.one.list))
	FIM[idx.zeta, idx.gamma] <- t(W) %*% (D * p*(1-p) * nu * S)
	FIM[idx.gamma, idx.zeta] <- t(FIM[idx.zeta, idx.gamma])

	# FIM[beta, zeta]
	D <- unlist(Map(function(x) { x[1,3] }, FIM.one.list))
	FIM[idx.zeta, idx.beta] <- t(W) %*% (D * lambda * p*(1-p) * X)
	FIM[idx.beta, idx.zeta] <- t(FIM[idx.zeta, idx.beta])

	return(FIM)
}

fim.zicmp.old <- function(lambda, nu, p, max = 100)
{
	z <- z_hybrid(lambda, nu)
	dzdlambda <- z_prodj(lambda, nu, max)/lambda
	dzdnu <- -z_prodlogj(lambda, nu, max)
	d2zdlambda2 <- (1/lambda^2)*(z_prodj2(lambda, nu, max) -
		z_prodj(lambda, nu, max))
	d2zdnu2 <- z_prodlogj2(lambda, nu, max)
	d2zdlambdadnu <- -(1/lambda)*z_prodjlogj(lambda, nu, max)

	dlogzdlambda <- dzdlambda/z
	d2logzdlambda2 <- d2zdlambda2/z - dzdlambda^2/z^2
	dlogzdnu <- dzdnu / z
	d2logzdnu2 <- d2zdnu2/z - dzdnu^2 / z^2
	d2logzdlambdadnu <- d2zdlambdadnu / z - dzdnu * dzdlambda / z^2
	meany <- (1-p)*lambda*dlogzdlambda

	FIM <- matrix(NA, 3, 3)
	colnames(FIM) <- c("lambda", "nu", "p")
	rownames(FIM) <- c("lambda", "nu", "p")

	# FIM[lambda, lambda]
	FIM[1,1] <- (1-p)*d2logzdlambda2 - p*(1-p)*dlogzdlambda^2 / (p*(z-1)+1) + meany / lambda^2

	# FIM[nu, nu]
	FIM[2,2] <- (1-p)*d2logzdnu2 - p*(1-p)*dlogzdnu^2 / (p*(z-1)+1)
	
	# FIM[p, p]
	FIM[3,3] <- (1/z) * (z-1)^2 / (p*(z-1) + 1) + (1 - 1/z) / (1-p)

	# FIM[lambda, nu]
	FIM[1,2] <- (1-p)*d2logzdlambdadnu - p*(1-p)*dlogzdnu*dlogzdlambda / (p*(z-1)+1)
	FIM[2,1] <- FIM[1,2]

	# FIM[lambda, p]
	FIM[1,3] <- -1/(p*(z-1) + 1) * dlogzdlambda
	FIM[3,1] <- FIM[1,3]

	# FIM[nu, p]
	FIM[2,3] <- -1/(p*(z-1) + 1) * dlogzdnu
	FIM[3,2] <- FIM[2,3]

	return(FIM)
}

fim.zicmp.reg.old <- function(X, S, W, beta, gamma, zeta, max = 100)
{
	n <- nrow(X)
	qq <- ncol(X) + ncol(S) + ncol(W)
	lambda <- as.numeric(exp(X %*% beta))
	nu <- as.numeric(exp(S %*% gamma))
	p <- as.numeric(plogis(W %*% zeta))

	z <- z_hybrid(lambda, nu)
	dzdlambda <- z_prodj(lambda, nu, max) / lambda
	dzdnu <- -z_prodlogj(lambda, nu, max)
	d2zdlambda2 <- (1/lambda^2)*(z_prodj2(lambda, nu, max) -
		z_prodj(lambda, nu, max))
	d2zdnu2 <- z_prodlogj2(lambda, nu, max)
	d2zdlambdadnu <- -(1/lambda)*z_prodjlogj(lambda, nu, max)

	dlogzdlambda <- dzdlambda/z
	d2logzdlambda2 <- d2zdlambda2/z - dzdlambda^2/z^2
	dlogzdnu <- dzdnu / z
	d2logzdnu2 <- d2zdnu2/z - dzdnu^2 / z^2
	d2logzdlambdadnu <- d2zdlambdadnu / z - dzdnu * dzdlambda / z^2
	meany <- (1-p)*lambda*dlogzdlambda

	FIM <- matrix(NA, qq, qq)
	colnames(FIM) <- c(sprintf("beta%d", 1:length(beta)),
		sprintf("gamma%d", 1:length(gamma)),
		sprintf("zeta%d", 1:length(zeta)))
	rownames(FIM) <- colnames(FIM)

	idx.beta <- 1:length(beta)
	idx.gamma <- 1:length(gamma) + length(beta)
	idx.zeta <- 1:length(zeta) + length(gamma) + length(beta)

	# FIM[beta, beta]
	D <- (1-p)*d2logzdlambda2 - p*(1-p)*dlogzdlambda^2 / (p*(z-1)+1) + meany / lambda^2
	FIM[idx.beta, idx.beta] <- t(X) %*% ((D * lambda^2) * X)

	# FIM[gamma, gamma]
	D <- (1-p)*d2logzdnu2 - p*(1-p)*dlogzdnu^2 / (p*(z-1)+1)
	FIM[idx.gamma, idx.gamma] <- t(S) %*% ((D * nu^2) * S)

	# FIM[zeta, zeta]
	D <- (1/z) * (z-1)^2 / (p*(z-1) + 1) + (1 - 1/z) / (1-p)
	FIM[idx.zeta, idx.zeta] <- t(W) %*% ((D * p^2*(1-p)^2) * W)

	# FIM[beta, gamma]
	D <- (1-p)*d2logzdlambdadnu - p*(1-p)*dlogzdlambda*dlogzdnu / (p*(z-1) + 1)
	FIM[idx.beta, idx.gamma] <- t(X) %*% (D * lambda * nu * S)
	FIM[idx.gamma, idx.beta] <- t(FIM[idx.beta, idx.gamma])

	# FIM[zeta, gamma]
	D <- -1/(p*(z-1) + 1) * dlogzdnu
	FIM[idx.zeta, idx.gamma] <- t(W) %*% (D * p*(1-p) * nu * S)
	FIM[idx.gamma, idx.zeta] <- t(FIM[idx.zeta, idx.gamma])

	# FIM[beta, zeta]
	D <- -1/(p*(z-1) + 1) * dlogzdlambda
	FIM[idx.zeta, idx.beta] <- t(W) %*% (D * lambda * p*(1-p) * X)
	FIM[idx.beta, idx.zeta] <- t(FIM[idx.zeta, idx.beta])

	return(FIM)
}
