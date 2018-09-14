#include <Rcpp.h>
#include "cmp.h"
#include "util.h"

// Assume that quantiles q are not given on the log scale
// [[Rcpp::export]]
Rcpp::NumericVector qzicmp_cpp(const Rcpp::NumericVector& q,
	const Rcpp::NumericVector& lambda, const Rcpp::NumericVector& nu,
	const Rcpp::NumericVector& p, double tol = 1e-6)
{
	unsigned int n = q.size();
	if (n != lambda.size() || n != nu.size() || n != p.size())  {
		Rcpp::stop("lambda, nu, and p must have length n");
	}

	Rcpp::NumericVector x(n);
	for (unsigned int i = 0; i < n; i++) {
		Rcpp::NumericVector cmp_probs = cmp_allprobs(lambda(i), nu(i), tol);
		Rcpp::NumericVector zicmp_probs = (1 - p(i))*cmp_probs;
		zicmp_probs(0) += p(i);
		x(i) = qdiscrete(q(i), zicmp_probs);
	}

	return x;
}
