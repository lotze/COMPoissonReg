rqres.zicmp <- function(y, lambda, nu, p)
{
	n <- length(y)
	if (length(lambda) == 1) lambda <- rep(lambda, n)
	if (length(nu) == 1) nu <- rep(nu, n)
	if (length(p) == 1) p <- rep(p, n)

	F <- function(y) {
		ret <- numeric(n)
		for (i in 1:n) {
			ret[i] <- p.zi.compoisson(y[i], lambda[i], nu[i], p[i])
		}
		return(ret)
	}
	rqres(y, F)
}

rqres.cmp <- function(y, lambda, nu)
{
	n <- length(y)
	if (length(lambda) == 1) lambda <- rep(lambda, n)

	F <- function(y) {
		ret <- numeric(n)
		for (i in 1:n) {
			ret[i] <- pcom(y[i], lambda[i], nu)
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
