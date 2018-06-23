#include <Rcpp.h>
#include <list>
#include "COMPoissonReg.h"

// [[Rcpp::export]]
Rcpp::NumericVector dcmp_cpp(const Rcpp::NumericVector& x,
	const Rcpp::NumericVector& lambda, const Rcpp::NumericVector& nu,
	double tol = 1e-6, bool take_log = false)
{
	unsigned int n = x.size();
	if (n != lambda.size() || n != nu.size()) {
		Rcpp::stop("lambda and nu must have length n");
	}

	Rcpp::NumericVector fx(n);
	fx.fill(0);

	for (unsigned int i = 0; i < n; i++) {
		Rcpp::NumericVector allprobs = cmp_allprobs(lambda(i), nu(i), tol);

		// Rcpp::print(allprobs);
		// Rprintf("Looking for element with index %d", x(i) + 1);

		if (x(i) <= allprobs.size() - 1) {
			fx(i) = allprobs(x(i));
		}
	}

	if (take_log) {
		return log(fx);
	} else {
		return fx;
	}
}

// [[Rcpp::export]]
Rcpp::NumericVector pcmp_cpp(const Rcpp::NumericVector& x,
	const Rcpp::NumericVector& lambda, const Rcpp::NumericVector& nu,
	double tol = 1e-6)
{
	unsigned int n = x.size();
	if (n != lambda.size() || n != nu.size()) {
		Rcpp::stop("lambda and nu must have length n");
	}

	Rcpp::NumericVector Fx(n);
	Fx.fill(1);

	for (unsigned int i = 0; i < n; i++) {
		Rcpp::NumericVector allprobs = cmp_allprobs(lambda(i), nu(i), tol);
		Rcpp::NumericVector cumprobs = Rcpp::cumsum(allprobs);
		if (x(i) <= allprobs.size() - 1) {
			Fx(i) = cumprobs(x(i));
		}
	}

	return Fx;
}

// Assume that quantiles q are not given on the log scale
// [[Rcpp::export]]
Rcpp::NumericVector qcmp_cpp(const Rcpp::NumericVector& q,
	const Rcpp::NumericVector& lambda, const Rcpp::NumericVector& nu,
	double tol = 1e-6)
{
	unsigned int n = q.size();
	if (n != lambda.size() || n != nu.size()) {
		Rcpp::stop("lambda and nu must have length n");
	}

	Rcpp::NumericVector x(n);
	for (unsigned int i = 0; i < n; i++) {
		Rcpp::NumericVector allprobs = cmp_allprobs(lambda(i), nu(i), tol);
		x(i) = qdiscrete(q(i), allprobs);
	}

	return x;
}

// Produce n(i) iid draws for each lambda(i) and nu(i).
// If n is a scalar, take n(i) to be n
// [[Rcpp::export]]
Rcpp::NumericVector rcmp_cpp(unsigned int n, const Rcpp::NumericVector& lambda,
	const Rcpp::NumericVector& nu, double tol = 1e-6)
{
	if (n != lambda.size() || n != nu.size()) {
		Rcpp::stop("lambda and nu must have length n");
	}

	Rcpp::NumericVector u = Rcpp::runif(n, 0.0, 1.0);
	return qcmp_cpp(u, lambda, nu, tol);
}

// [[Rcpp::export]]
Rcpp::NumericVector cmp_expected_value(const Rcpp::NumericVector& lambda,
	const Rcpp::NumericVector& nu, double tol = 1e-6)
{
	unsigned int n = lambda.size();
	if (n != nu.size()) {
		Rcpp::stop("lambda and nu must be the same length");
	}
	
	Rcpp::NumericVector ex(n);
	ex.fill(0);

	for (unsigned int i = 0; i < n; i++) {
		Rcpp::NumericVector allprobs = cmp_allprobs(lambda(i), nu(i), tol);
		for (unsigned int j = 0; j < allprobs.size(); j++) {
			ex(i) += j * allprobs(j);
		}
	}

	return ex;
}
