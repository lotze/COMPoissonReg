dcmp <- function (x, lambda, nu, log = FALSE)
{
	prep <- prep.zicmp(length(x), lambda, nu)
	dcmp_cpp(x, prep$lambda, prep$nu, take_log = log)
}

rcmp <- function(n, lambda, nu)
{
	prep <- prep.zicmp(n, lambda, nu)
	rcmp_cpp(n, prep$lambda, prep$nu)
}

pcmp <- function(x, lambda, nu)
{
	prep <- prep.zicmp(length(x), lambda, nu)
	pcmp_cpp(x, prep$lambda, prep$nu)
}

qcmp <- function(q, lambda, nu, log.p = FALSE)
{
	prep <- prep.zicmp(length(q), lambda, nu)
	if (log.p) {
		log.q <- q
	} else {
		log.q <- log(q)
	}
	qcmp_cpp(log.q, prep$lambda, prep$nu)
}

cmp_expected_value <- function(lambda, nu)
{
	n <- max(length(lambda), length(nu))
	prep <- prep.zicmp(n, lambda, nu)

	res <- numeric(n)
	for (i in 1:n) {
		res[i] <- prep$lambda[i] * grad.fwd(z_hybrid, prep$lambda[i], nu = prep$nu[i], take_log = TRUE)
	}
	return(res)
}
