r.zi.compoisson <- function(n, lambda, nu, p)
{
	x <- integer(n)
	z <- rbinom(n, size = 1, prob = p)
	x[z == 1] <- 0
	x[z == 0] <- rcom(sum(z == 0), lambda, nu)
	return(x)
}

# Adapted from compoisson package
d.compoisson <- function (x, lambda, nu, log = FALSE) 
{
	if (lambda < 0 || nu < 0) 
		stop("Invalid arguments, only defined for lambda >= 0, nu >= 0")

	log.ff <- log(x >= 0 & (x == floor(x))) + x*log(lambda) - nu*lgamma(x+1) -
		com.compute.log.z(lambda, nu)
	if (log) log.ff
	else exp(log.ff)
}

# Warning: This should be improved for efficiency and convenience
p.compoisson <- function(x, lambda, nu)
{
	if (length(x) > 1)
		stop("Currently, x must be a single integer")
	if (length(length(lambda)) > 1)
		stop("Currently, lambda must be a single number")
	if (length(length(nu)) > 1)
		stop("Currently, nu must be a single number")

	if (x < 0)
		return(0)
	else
		return(sum(d.compoisson(0:floor(x), lambda, nu)))

}

d.zi.compoisson <- function(x, lambda, nu, p, log = FALSE)
{
	fx <- p*(x==0) + (1-p)*d.compoisson(x, lambda, nu)
	if (log) return(log(fx))
	else return(fx)
}

p.zi.compoisson <- function(x, lambda, nu, p)
{
	p*(x >= 0) + (1-p)*p.compoisson(x, lambda, nu)
}

