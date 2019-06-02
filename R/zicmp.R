dzicmp = function(x, lambda, nu, p, log = FALSE)
{
	prep = prep.zicmp(length(x), lambda, nu, p)
	fx = prep$p*(x==0) + (1-prep$p)*dcmp(x, prep$lambda, prep$nu)
	if (log) {
		return(log(fx))
	} else {
		return(fx)
	}
}

rzicmp = function(n, lambda, nu, p)
{
	prep = prep.zicmp(n, lambda, nu, p)
	s = rbinom(n, size = 1, prob = prep$p)
	(1-s) * rcmp(n, prep$lambda, prep$nu)
}

pzicmp = function(x, lambda, nu, p)
{
	prep = prep.zicmp(length(x), lambda, nu, p)
	prep$p*(x >= 0) + (1-prep$p)*pcmp(x, prep$lambda, prep$nu)
}

qzicmp = function(q, lambda, nu, p, log.p = FALSE)
{
	prep = prep.zicmp(length(q), lambda, nu, p)
	if (log.p) {
		log.q = q
	} else {
		log.q = log(q)
	}
	qzicmp_cpp(log.q, prep$lambda, prep$nu, prep$p)
}

zicmp_expected_value = function(lambda, nu, p)
{
	(1-p) * cmp_expected_value(lambda, nu)
}

# Extend lambda, nu, and p vectors to be compatible lengths.
# If all are length 1, do not extend them - this is a special case which
# is handled more efficiently.
prep.zicmp = function(n, lambda, nu, p = 0)
{
	L = max(length(lambda), length(nu), length(p))

	if (n > 1 && L > 1) { stopifnot(n == L) }
	if (length(lambda) == 1 && L > 1) { lambda = rep(lambda, L) }
	if (length(nu) == 1 && L > 1) { nu = rep(nu, L) }
	if (length(p) == 1 && L > 1) { p = rep(p, L) }

	list(lambda = lambda, nu = nu, p = p)
}
