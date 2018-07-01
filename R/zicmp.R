rzicmp <- function(n, lambda, nu, p)
{
	if (length(lambda) == 1) { lambda <- rep(lambda,n) }
	if (length(nu) == 1) { nu <- rep(nu,n) }
	if (length(p) == 1) { p <- rep(p,n) }
	x <- integer(n)
	s <- rbinom(n, size = 1, prob = p)
	idx <- which(s == 0)
	x[idx] <- rcmp(length(idx), lambda[idx], nu[idx])
	return(x)
}

dzicmp <- function(x, lambda, nu, p, log = FALSE)
{
	n <- length(x)
	if (length(lambda) == 1) { lambda <- rep(lambda,n) }
	if (length(nu) == 1) { nu <- rep(nu,n) }
	if (length(p) == 1) { p <- rep(p,n) }
	fx <- p*(x==0) + (1-p)*dcmp(x, lambda, nu)
	if (log) return(log(fx))
	else return(fx)
}

pzicmp <- function(x, lambda, nu, p)
{
	n <- length(x)
	if (length(lambda) == 1) { lambda <- rep(lambda,n) }
	if (length(nu) == 1) { nu <- rep(nu,n) }
	if (length(p) == 1) { p <- rep(p,n) }
	p*(x >= 0) + (1-p)*pcmp(x, lambda, nu)
}

qzicmp <- function(q, lambda, nu, p, log.p = FALSE)
{
	n <- length(q)
	log.q <- ifelse(log.p, q, log(q))
	if (length(lambda) == 1) { lambda <- rep(lambda, n) }
	if (length(nu) == 1) { nu <- rep(nu, n) }
	if (length(p) == 1) { p <- rep(p, n) }
	qzicmp_cpp(q, lambda, nu, p)
}

zicmp_expected_value <- function(lambda, nu, p)
{
	(1-p) * cmp_expected_value(lambda, nu)
}
