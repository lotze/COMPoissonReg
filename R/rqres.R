rqres.zicmp <- function(y, lambda, nu, p)
{
	n <- length(y)
	if (length(lambda) == 1) lambda <- rep(lambda, n)
	if (length(p) == 1) p <- rep(p, n)

	F <- function(y) {
		ret <- numeric(n)
		for (i in 1:n) {
			ret[i] <- p.zi.compoisson(y[i], lambda[i], nu, p[i])
		}
		return(ret)
	}
	rqres(y, F)
}

rqres.zicmp.reg <- function(y, X, W, beta, nu, zeta)
{
	n <- length(y)
	lambda <- exp(X %*% beta)
	p <- plogis(W %*% zeta)

	F <- function(y) {
		ret <- numeric(n)
		for (i in 1:n) {
			ret[i] <- p.zi.compoisson(y[i], lambda[i], nu, p[i])
		}
		return(ret)
	}

	rqres(y, F)
}

rqres <- function(y, F, eps = 1e-6)
{
	n <- length(y)
	FL <- F(y - eps)
	FU <- F(y)
	u <- runif(n, min = FL, max = FU)
	qres <- qnorm(u)
	return(qres)
}
