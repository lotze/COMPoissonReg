# Warning: This should be improved for efficiency and convenience
pcom <- function(x, lambda, nu, max = 100)
{
	if (length(x) > 1)
		stop("Currently, x must be a single integer")
	if (length(length(lambda)) > 1)
		stop("Currently, lambda must be a single number")
	if (length(length(nu)) > 1)
		stop("Currently, nu must be a single number")

	if (x < 0) {
		return(0)
	}
	else {
		z <- computez(lambda, nu, max)
		return(sum(dcom(0:floor(x), lambda, nu, z)))
	}
}

r.zi.compoisson <- function(n, lambda, nu, p)
{
	if (length(lambda) == 1) { lambda <- rep(lambda,n) }
	if (length(nu) == 1) { nu <- rep(nu,n) }
	if (length(p) == 1) { p <- rep(p,n) }

	x <- integer(n)
	s <- rbinom(n, size = 1, prob = p)
	x[s == 1] <- 0

	# rcom gives warnings if we try to do a vectorized call.
	# Just use a loop for now.
	for (i in which(s == 0)) {
		x[i] <- rcom(1, lambda[i], nu[i])
	}

	return(x)
}

d.zi.compoisson <- function(x, lambda, nu, p, max = 100, log = FALSE)
{
	z <- computez(lambda, nu, max)
	fx <- p*(x==0) + (1-p)*dcom(x, lambda, nu, z)
	if (log) return(log(fx))
	else return(fx)
}

p.zi.compoisson <- function(x, lambda, nu, p, max = 100)
{
	p*(x >= 0) + (1-p)*pcom(x, lambda, nu, max = max)
}
