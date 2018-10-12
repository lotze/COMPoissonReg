rqres.zicmp <- function(y, lambda, nu, p)
{
	F <- function(y) {
		pzicmp(y, lambda, nu, p)
	}
	rqres(y, F)
}

rqres.cmp <- function(y, lambda, nu)
{
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
