dzip = function(x, lambda, p, log = FALSE)
{
	fx = p*(x==0) + (1-p)*dpois(x, lambda)
	if (log) return(log(fx))
	else return(fx)
}
