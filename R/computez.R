# Sum (from j=0 to j=max) of lambda^j/((j!)^nu)
computez <- function(lambda, nu, max)
{
	n <- length(lambda)
	L <- matrix(log(lambda), nrow=n, ncol=max+1, byrow = FALSE)
	M <- matrix(nu, nrow=n, ncol=max+1, byrow = FALSE)
	J <- matrix(0:max, nrow=n, ncol=max+1, byrow = TRUE)
	log.res <- J*L - M*lgamma(J+1)
	rowSums(exp(log.res))
}

# Sum (from j=0 to j=max) of j*lambda^j/((j!)^nu)
computez.prodj <- function(lambda, nu, max)
{
	n <- length(lambda)
	L <- matrix(log(lambda), nrow=n, ncol=max+1, byrow = FALSE)
	M <- matrix(nu, nrow=n, ncol=max+1, byrow = FALSE)
	J <- matrix(0:max, nrow=n, ncol=max+1, byrow = TRUE)
	log.res <- log(J) + J*L - M*lgamma(J+1)
	rowSums(exp(log.res))
}

# Sum (from j=0 to j=max) of j^2*lambda^j/((j!)^nu)
computez.prodj2 <- function(lambda, nu, max)
{
	n <- length(lambda)
	L <- matrix(log(lambda), nrow=n, ncol=max+1, byrow = FALSE)
	M <- matrix(nu, nrow=n, ncol=max+1, byrow = FALSE)
	J <- matrix(0:max, nrow=n, ncol=max+1, byrow = TRUE)
	log.res <- 2*log(J) + J*L - M*lgamma(J+1)
	rowSums(exp(log.res))
}

# Sum (from j=0 to j=max) of jlog(j!)*lambda^j/((j!)^nu)
computez.prodjlogj <- function(lambda, nu, max)
{
	n <- length(lambda)
	L <- matrix(log(lambda), nrow=n, ncol=max+1, byrow = FALSE)
	M <- matrix(nu, nrow=n, ncol=max+1, byrow = FALSE)
	J <- matrix(0:max, nrow=n, ncol=max+1, byrow = TRUE)
	log.res <- log(J) + log(lgamma(J+1)) + J*L - M*lgamma(J+1)
	rowSums(exp(log.res))
}

# Sum (from j=0 to j=max) of log(j!)*lambda^j/((j!)^nu)
computez.prodlogj <- function(lambda, nu, max)
{
	n <- length(lambda)
	L <- matrix(log(lambda), nrow=n, ncol=max+1, byrow = FALSE)
	M <- matrix(nu, nrow=n, ncol=max+1, byrow = FALSE)
	J <- matrix(0:max, nrow=n, ncol=max+1, byrow = TRUE)
	log.res <- log(lgamma(J+1)) + J*L - M*lgamma(J+1)
	rowSums(exp(log.res))
}

# Sum (from j=0 to j=max) of (log(j!))^2 * lambda^j/((j!)^nu)
computez.prodlogj2 <- function(lambda, nu, max)
{
	n <- length(lambda)
	L <- matrix(log(lambda), nrow=n, ncol=max+1, byrow = FALSE)
	M <- matrix(nu, nrow=n, ncol=max+1, byrow = FALSE)
	J <- matrix(0:max, nrow=n, ncol=max+1, byrow = TRUE)
	log.res <- 2*log(lgamma(J+1)) + J*L - M*lgamma(J+1)
	rowSums(exp(log.res))
}

