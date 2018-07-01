dzip <- function(x, lambda, p, log = FALSE)
{
	n <- length(x)
	if (length(lambda) == 1) { lambda <- rep(lambda,n) }
	if (length(p) == 1) { p <- rep(p,n) }
	fx <- p*(x==0) + (1-p)*dpois(x, lambda)
	if (log) return(log(fx))
	else return(fx)
}
