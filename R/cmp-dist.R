dcmp <- function (x, lambda, nu, z = NULL, log = FALSE, max = 100) 
{
	if (is.null(z)) {
		z <- computez(lambda, nu, max)
	}
	
	log.ff <- x * log(lambda) - nu * lgamma(x+1) - log(z)
	if (log) {
		return(log.ff)
	} else {
		return(exp(log.ff))
	}
}

rcmp <- function(n, lambda, nu, max = 100)
{
	if (length(lambda) == 1) { lambda <- rep(lambda, n) }
	if (length(nu) == 1) { nu <- rep(nu, n) }
	
	u <- runif(n)
	y <- numeric(n)
	z <- computez(lambda, nu, max)
	
	for (i in 1:n) {
		py <- 1 / z[i]
		while (py < u[i]) {
			y[i] <- y[i] + 1
			py <- py + dcmp(y[i], lambda[i], nu[i], z = z[i])
		}
	}

	return(y)
}

pcmp <- function(x, lambda, nu, max = 100, z = NULL)
{
	n <- length(x)
	if (length(lambda) == 1) { lambda <- rep(lambda, n) }
	if (length(nu) == 1) { nu <- rep(nu, n) }

	if (is.null(z)) {
		z <- computez(lambda, nu, max)
	}
	Fx <- numeric(n)
	
	for (i in 1:n) {
		if (x[i] > 0) {
			Fx[i] <- sum(dcmp(0:floor(x[i]), lambda[i], nu[i], z[i]))
		}
	}

	return(Fx)
}

rzicmp <- function(n, lambda, nu, p, max = 100)
{
	if (length(lambda) == 1) { lambda <- rep(lambda,n) }
	if (length(nu) == 1) { nu <- rep(nu,n) }
	if (length(p) == 1) { p <- rep(p,n) }

	x <- integer(n)
	s <- rbinom(n, size = 1, prob = p)
	x[s == 0] <- rcmp(sum(s == 0), lambda[s == 0], nu[s == 0], max = max)
	return(x)
}

dzicmp <- function(x, lambda, nu, p, max = 100, log = FALSE)
{
	z <- computez(lambda, nu, max)
	fx <- p*(x==0) + (1-p)*dcmp(x, lambda, nu, z)
	if (log) return(log(fx))
	else return(fx)
}

pzicmp <- function(x, lambda, nu, p, max = 100)
{
	p*(x >= 0) + (1-p)*pcmp(x, lambda, nu, max = max)
}
