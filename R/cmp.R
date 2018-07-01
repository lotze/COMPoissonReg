dcmp <- function (x, lambda, nu, log = FALSE)
{
	n <- length(x)
	if (length(lambda) == 1) { lambda <- rep(lambda, n) }
	if (length(nu) == 1) { nu <- rep(nu, n) }
	dcmp_cpp(x, lambda, nu, take_log = log)
}

rcmp <- function(n, lambda, nu)
{
	if (length(lambda) == 1) { lambda <- rep(lambda, n) }
	if (length(nu) == 1) { nu <- rep(nu, n) }
	rcmp_cpp(n, lambda, nu)
}

pcmp <- function(x, lambda, nu)
{
	n <- length(x)
	if (length(lambda) == 1) { lambda <- rep(lambda, n) }
	if (length(nu) == 1) { nu <- rep(nu, n) }
	pcmp_cpp(x, lambda, nu)
}

qcmp <- function(q, lambda, nu, log.p = FALSE)
{
	n <- length(q)
	log.q <- ifelse(log.p, q, log(q))
	if (length(lambda) == 1) { lambda <- rep(lambda, n) }
	if (length(nu) == 1) { nu <- rep(nu, n) }
	qcmp_cpp(x, lambda, nu)
}

cmp_expected_value <- function(lambda, nu)
{
	n <- max(length(lambda), length(nu))
	if (length(lambda) == 1) { lambda <- rep(lambda, n) }
	if (length(nu) == 1) { nu <- rep(nu, n) }
	res <- numeric(n)
	for (i in 1:n) {
		res[i] <- lambda[i] * grad.fwd(z_hybrid, lambda[i], nu = nu[i], take_log = TRUE)
	}
	return(res)
}
