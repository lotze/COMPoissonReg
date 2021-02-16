#include "z-derivs.h"
#include "util.h"

Rcpp::NumericVector z_prodj(const Rcpp::NumericVector& lambda,
	const Rcpp::NumericVector& nu, unsigned int max)
{
	unsigned int n = lambda.size();
	if (n != nu.size()) {
		Rcpp::stop("lambda and nu must be the same length");
	}

	Rcpp::NumericVector res(n);
	res.fill(0);

	for (unsigned int j = 0; j < max; j++) {
		res += exp(log(j) + j*log(lambda) - nu*lgamma(j+1));
	}

	return res;
}

Rcpp::NumericVector z_prodj2(const Rcpp::NumericVector& lambda,
	const Rcpp::NumericVector& nu, unsigned int max)
{
	unsigned int n = lambda.size();
	if (n != nu.size()) {
		Rcpp::stop("lambda and nu must be the same length");
	}

	Rcpp::NumericVector res(n);
	res.fill(0);

	for (unsigned int j = 0; j < max; j++) {
		res += exp(2*log(j) + j*log(lambda) - nu*lgamma(j+1));
	}

	return res;
}

Rcpp::NumericVector z_prodjlogj(const Rcpp::NumericVector& lambda,
	const Rcpp::NumericVector& nu, unsigned int max)
{
	unsigned int n = lambda.size();
	if (n != nu.size()) {
		Rcpp::stop("lambda and nu must be the same length");
	}

	Rcpp::NumericVector res(n);
	res.fill(0);

	for (unsigned int j = 0; j < max; j++) {
		res += exp(log(j) + log(lgamma(j+1)) + j*log(lambda) - nu*lgamma(j+1));
	}

	return res;
}

Rcpp::NumericVector z_prodlogj(const Rcpp::NumericVector& lambda,
	const Rcpp::NumericVector& nu, unsigned int max)
{
	unsigned int n = lambda.size();
	if (n != nu.size()) {
		Rcpp::stop("lambda and nu must be the same length");
	}

	Rcpp::NumericVector res(n);
	res.fill(0);

	for (unsigned int j = 0; j < max; j++) {
		res += exp(log(lgamma(j+1)) + j*log(lambda) - nu*lgamma(j+1));
	}

	return res;
}

Rcpp::NumericVector z_prodlogj2(const Rcpp::NumericVector& lambda,
	const Rcpp::NumericVector& nu, unsigned int max)
{
	unsigned int n = lambda.size();
	if (n != nu.size()) {
		Rcpp::stop("lambda and nu must be the same length");
	}

	Rcpp::NumericVector res(n);
	res.fill(0);

	for (unsigned int j = 0; j < max; j++) {
		res += exp(2*log(lgamma(j+1)) + j*log(lambda) - nu*lgamma(j+1));
	}

	return res;
}
