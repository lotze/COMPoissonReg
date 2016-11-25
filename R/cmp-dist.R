# Warning: This should be improved for efficiency and convenience
pcom <- function(x, lambda, nu)
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
		return(sum(dcom(0:floor(x), lambda, nu)))
}

r.zi.compoisson <- function(n, lambda, nu, p)
{
	x <- integer(n)
	z <- rbinom(n, size = 1, prob = p)
	x[z == 1] <- 0
	x[z == 0] <- rcom(sum(z == 0), lambda, nu)
	return(x)
}

d.zi.compoisson <- function(x, lambda, nu, p, log = FALSE)
{
	fx <- p*(x==0) + (1-p)*dcom(x, lambda, nu)
	if (log) return(log(fx))
	else return(fx)
}

p.zi.compoisson <- function(x, lambda, nu, p)
{
	p*(x >= 0) + (1-p)*pcom(x, lambda, nu)
}
