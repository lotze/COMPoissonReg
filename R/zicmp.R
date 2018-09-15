rzicmp <- function(n, lambda, nu, p)
{
	s <- rbinom(n, size = 1, prob = p)
	(1-s) * rcmp(n, lambda, nu)
}

dzicmp <- function(x, lambda, nu, p, log = FALSE)
{
	fx <- p*(x==0) + (1-p)*dcmp(x, lambda, nu)
	if (log) {
		return(log(fx))
	} else {
		return(fx)
	}
}

pzicmp <- function(x, lambda, nu, p)
{
	p*(x >= 0) + (1-p)*pcmp(x, lambda, nu)
}

qzicmp <- function(q, lambda, nu, p, log.p = FALSE)
{
	if (log.p) {
		log.q <- q
	} else {
		log.q <- log(q)
	}
	qzicmp_cpp(log.q, lambda, nu, p)
}

zicmp_expected_value <- function(lambda, nu, p)
{
	(1-p) * cmp_expected_value(lambda, nu)
}
