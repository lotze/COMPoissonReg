rqres.zicmp <- function(y, lambda, nu, p)
{
	n <- length(y)
	if (length(lambda) == 1) lambda <- rep(lambda, n)
	if (length(nu) == 1) nu <- rep(nu, n)
	if (length(p) == 1) p <- rep(p, n)

	F <- function(y) {
		pzicmp(y, lambda, nu, p)
	}
	rqres(y, F)
}

rqres.cmp <- function(y, lambda, nu)
{
	n <- length(y)
	if (length(lambda) == 1) lambda <- rep(lambda, n)

	F <- function(y) {
		pcmp(y, lambda, nu)
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
