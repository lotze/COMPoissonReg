# given estimates for beta, zeta, and nu, this code computes components of information matrix
# constant dispersion parameter, nu
# ln(lambda) = x %*% beta
# logit(p) = w %*% zeta

fim.zicmp <- function(lambda, nu, p)
{
	z <- z_hybrid(lambda, nu)
	dzdlambda <- computez.prodj(lambda, nu)/lambda
	dzdnu <- -computez.prodlogj(lambda, nu)
	d2zdlambda2 <- (1/lambda^2)*(computez.prodj2(lambda, nu) -
		computez.prodj(lambda, nu))
	d2zdnu2 <- computez.prodlogj2(lambda, nu)
	d2zdlambdadnu <- -(1/lambda)*computez.prodjlogj(lambda, nu)

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

fim.zicmp.reg <- function(X, S, W, beta, gamma, zeta)
{
	n <- nrow(X)
	qq <- ncol(X) + ncol(S) + ncol(W)
	lambda <- as.numeric(exp(X %*% beta))
	nu <- as.numeric(exp(S %*% gamma))
	p <- as.numeric(plogis(W %*% zeta))

	z <- computez(lambda, nu)
	dzdlambda <- computez.prodj(lambda, nu) / lambda
	dzdnu <- -computez.prodlogj(lambda, nu)
	d2zdlambda2 <- (1/lambda^2)*(computez.prodj2(lambda, nu) -
		computez.prodj(lambda, nu))
	d2zdnu2 <- computez.prodlogj2(lambda, nu)
	d2zdlambdadnu <- -(1/lambda)*computez.prodjlogj(lambda, nu)

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

var.zicmp.reg <- function(X, S, W, beta, gamma, zeta)
{
	FIM <- fim.zicmp.reg(X, S, W, beta, gamma, zeta)
	solve(FIM)
}
