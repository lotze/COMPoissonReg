# given estimates for beta, zeta, and nu, this code computes components of information matrix
# constant dispersion parameter, nu
# ln(lambda) = x %*% beta
# logit(p) = w %*% zeta

fim.zicmp <- function(lambda, nu, p, max = 100)
{
	z <- computez.lambdaest(lambda,nu,max=max)
	dzdlambda <- computez.prodj.lambdaest(lambda, nu, max=max)/lambda
	dzdnu <- -computez.prodlogj.lambdaest(lambda, nu, max=max)
	d2zdlambda2 <- (1/lambda^2)*(computez.prodj2.lambdaest(lambda, nu, max=max) -
		computez.prodj.lambdaest(lambda, nu, max=max))
	d2zdnu2 <- computez.prodlogj2.lambdaest(lambda, nu, max=max)
	d2zdlambdadnu <- -(1/lambda)*computez.prodjlogj.lambdaest(lambda, nu, max=max)

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

fim.zicmp.reg <- function(x, beta, w, zeta, nu, max = 100)
{
	n <- nrow(x)
	qq <- ncol(x) + 1 + ncol(w)
	lambda <- as.numeric(exp(x %*% beta))
	p <- as.numeric(plogis(w %*% zeta))

	z <- computez.lambdaest(lambda, nu, max=max)
	dzdlambda <- computez.prodj.lambdaest(lambda, nu, max=max) / lambda
	dzdnu <- -computez.prodlogj.lambdaest(lambda, nu, max=max)
	d2zdlambda2 <- (1/lambda^2)*(computez.prodj2.lambdaest(lambda, nu, max=max) -
		computez.prodj.lambdaest(lambda, nu, max=max))
	d2zdnu2 <- computez.prodlogj2.lambdaest(lambda, nu, max=max)
	d2zdlambdadnu <- -(1/lambda)*computez.prodjlogj.lambdaest(lambda, nu, max=max)

	dlogzdlambda <- dzdlambda/z
	d2logzdlambda2 <- d2zdlambda2/z - dzdlambda^2/z^2
	dlogzdnu <- dzdnu / z
	d2logzdnu2 <- d2zdnu2/z - dzdnu^2 / z^2
	d2logzdlambdadnu <- d2zdlambdadnu / z - dzdnu * dzdlambda / z^2
	meany <- (1-p)*lambda*dlogzdlambda

	FIM <- matrix(NA, qq, qq)
	colnames(FIM) <- c(sprintf("beta%d", 1:length(beta)), "nu", sprintf("zeta%d", 1:length(zeta)))
	rownames(FIM) <- colnames(FIM)

	idx.beta <- 1:length(beta)
	idx.nu <- 1 + length(beta)
	idx.zeta <- 1:length(zeta) + 1 + length(beta)

	# FIM[beta, beta]
	D <- (1-p)*d2logzdlambda2 - p*(1-p)*dlogzdlambda^2 / (p*(z-1)+1) + meany / lambda^2
	FIM[idx.beta, idx.beta] <- t(x) %*% ((D * lambda^2) * x)

	# FIM[nu, nu]
	D <- (1-p)*d2logzdnu2 - p*(1-p)*dlogzdnu^2 / (p*(z-1)+1)
	FIM[idx.nu, idx.nu] <- sum(D)

	# FIM[zeta, zeta]
	D <- (1/z) * (z-1)^2 / (p*(z-1) + 1) + (1 - 1/z) / (1-p)
	FIM[idx.zeta, idx.zeta] <- t(w) %*% ((D * p^2*(1-p)^2) * w)

	# FIM[beta, nu]
	D <- (1-p)*d2logzdlambdadnu - p*(1-p)*dlogzdlambda*dlogzdnu / (p*(z-1) + 1)
	FIM[idx.beta, idx.nu] <- matrix(1, 1, n) %*% ((D * lambda) * x)
	FIM[idx.nu, idx.beta] <- FIM[idx.beta, idx.nu]

	# FIM[zeta, nu]
	D <- -1/(p*(z-1) + 1) * dlogzdnu
	FIM[idx.zeta, idx.nu] <- matrix(1, 1, n) %*% (D * p*(1-p) * w)
	FIM[idx.nu, idx.zeta] <- FIM[idx.zeta, idx.nu]

	# FIM[beta, zeta]
	D <- -1/(p*(z-1) + 1) * dlogzdlambda
	FIM[idx.zeta, idx.beta] <- t(w) %*% (D * lambda * p*(1-p) * x)
	FIM[idx.beta, idx.zeta] <- t(FIM[idx.zeta, idx.beta])
	
	return(FIM)
}

var.zicmp.reg <- function(x, beta, w, zeta, nu, max = 100)
{
	FIM <- fim.zicmp.reg(x, beta, w, zeta, nu, max = max)
	solve(FIM)
}
