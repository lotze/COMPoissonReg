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
	x <- numeric(n)
	logz <- computez(lambda, nu, max, log = TRUE, autoscale = TRUE)

	for (i in 1:n) {
		px <- exp(x[i]*log(lambda[i]) - nu[i]*lgamma(x[i]+1) - logz[i])
			# dcmp(x[i], lambda[i], nu[i], z = z[i], max = max)
		while (px < u[i]) {
			x[i] <- x[i] + 1
			px <- px + exp(x[i]*log(lambda[i]) - nu[i]*lgamma(x[i]+1) - logz[i])
				# dcmp(x[i], lambda[i], nu[i], z = z[i], max = max)
			print(px)
		}
	}

	return(x)
}

pcmp <- function(x, lambda, nu, max = 100)
{
	n <- length(x)
	if (length(lambda) == 1) { lambda <- rep(lambda, n) }
	if (length(nu) == 1) { nu <- rep(nu, n) }

	z <- computez(lambda, nu, max)
	Fx <- numeric(n)

	for (i in 1:n) {
		if (x[i] > 0) {
			Fx[i] <- sum(dcmp(0:floor(x[i]), lambda[i], nu[i], z = z[i], max = max))
		}
	}

	return(Fx)
}

qcmp <- function(q, lambda, nu, max = 100, log.p = FALSE)
{
	n <- length(q)
	if (length(log.p) == 1) { log.p <- rep(log.p, n) }
	log.q <- ifelse(log.p, q, log(q))
	if (length(lambda) == 1) { lambda <- rep(lambda, n) }
	if (length(nu) == 1) { nu <- rep(nu, n) }

	x <- numeric(n)
	z <- computez(lambda, nu, max)

	for (i in 1:n) {
		px <- dcmp(x[i], lambda = lambda[i], nu = nu[i], z = z[i], max = max)
		while (log(px) < log.q[i]) {
			x[i] <- x[i] + 1
			px <- px + dcmp(x[i], lambda = lambda[i], nu = nu[i], z = z[i], max = max)
		}
	}

	return(x)
}

rzicmp <- function(n, lambda, nu, p, max = 100)
{
	if (length(lambda) == 1) { lambda <- rep(lambda,n) }
	if (length(nu) == 1) { nu <- rep(nu,n) }
	if (length(p) == 1) { p <- rep(p,n) }

	x <- integer(n)
	s <- rbinom(n, size = 1, prob = p)
	idx <- which(s == 0)
	x[idx] <- rcmp(length(idx), lambda[idx], nu[idx], max = max)
	return(x)
}

dzicmp <- function(x, lambda, nu, p, z = NULL, max = 100, log = FALSE)
{
	if (is.null(z)) {
		z <- computez(lambda, nu, max)
	}
	fx <- p*(x==0) + (1-p)*dcmp(x, lambda, nu, z = z, max = max)
	if (log) return(log(fx))
	else return(fx)
}

pzicmp <- function(x, lambda, nu, p, max = 100)
{
	p*(x >= 0) + (1-p)*pcmp(x, lambda, nu, max = max)
}

qzicmp <- function(q, lambda, nu, p, max = 100, log.p = FALSE)
{
	n <- length(q)
	log.q <- ifelse(log.p, q, log(q))
	if (length(lambda) == 1) { lambda <- rep(lambda, n) }
	if (length(nu) == 1) { nu <- rep(nu, n) }
	if (length(p) == 1) { p <- rep(p, n) }

	x <- numeric(n)
	z <- computez(lambda, nu, max)

	for (i in 1:n) {
		px <- dzicmp(x[i], lambda = lambda[i], nu = nu[i], p = p[i], z = z[i], max = max)
		while (log(px) < log.q[i]) {
			x[i] <- x[i] + 1
			px <- px + dzicmp(x[i], lambda = lambda[i], nu = nu[i], p = p[i], z = z[i], max = max)
		}
	}
	
	return(x)
}
