#' COM-Poisson Distribution
#' 
#' Functions for the COM-Poisson distribution.
#' 
#' @param x vector of quantiles.
#' @param q vector of probabilities.
#' @param n number of observations.
#' @param lambda rate parameter.
#' @param p zero-inflation probability parameter.
#' @param log logical; if TRUE, probabilities are returned on log-scale.
#' @param log.p logical; if TRUE, probabilities \code{p} are given as \eqn{\log(p)}.
#' 
#' @return 
#' \describe{
#' \item{dzip}{the density,}
#' \item{pzip}{the cumulative probability,}
#' \item{qzip}{the quantile function,}
#' \item{rzip}{generates random values,}
#' \item{ezip}{the expected value,}
#' \item{vzip}{the variance,}
#' }
#' 
#' @author Kimberly Sellers
#' @keywords Zero-Inflated Poisson distribution
#' @name ZIP Distribution
NULL

#' @name ZIP Distribution
#' @export
dzip = function(x, lambda, p, log = FALSE)
{
	fx = p*(x==0) + (1-p)*dpois(x, lambda)
	if (log) { return(log(fx)) } else { return(fx) }
}

#' @name ZIP Distribution
#' @export
rzip = function(n, lambda, p)
{
	s = rbinom(n, size = 1, prob = p)
	(1-s) * rpois(n, lambda)
}

#' @name ZIP Distribution
#' @export
pzip = function(x, lambda, p)
{
	p*(x >= 0) + (1-p)*ppois(x, lambda)
}

#' @name ZIP Distribution
#' @export
qzip = function(q, lambda, p, log.p = FALSE)
{
	q_tx = (q-p) / (1-p)
	qpois(q_tx, lambda, log.p = log.p)
}

#' @name ZIP Distribution
#' @export
ezip = function(lambda, p)
{
	(1-p) * lambda
}

#' @name ZIP Distribution
#' @export
vzip = function(lambda, p)
{
	epois = lambda
	vpois = lambda
	(1-p) * (vpois + p*epois^2)
}
