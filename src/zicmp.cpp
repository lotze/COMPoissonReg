#include <Rcpp.h>
#include "cmp.h"
#include "util.h"

// Assume that quantiles q are given on the log scale
// [[Rcpp::export]]
Rcpp::NumericVector qzicmp_cpp(const Rcpp::NumericVector& logq,
	const Rcpp::NumericVector& lambda, const Rcpp::NumericVector& nu,
	const Rcpp::NumericVector& p, double tol = 1e-6)
{
	unsigned int n = logq.size();
	Rcpp::NumericVector x(n);
	
	if (lambda.size() == 1 && nu.size() == 1 &&  p.size() == 1)  {
		Rcpp::NumericVector cmp_probs = cmp_allprobs(lambda(0), nu(0), tol);
		Rcpp::NumericVector zicmp_probs = (1 - p(0))*cmp_probs;
		zicmp_probs(0) += p(0);
		for (unsigned int i = 0; i < n; i++) {
			x(i) = qdiscrete(exp(logq(i)), zicmp_probs);
		}
	} else if (lambda.size() == n && nu.size() == n && p.size() == n) {
		for (unsigned int i = 0; i < n; i++) {
			Rcpp::NumericVector cmp_probs = cmp_allprobs(lambda(i), nu(i), tol);
			Rcpp::NumericVector zicmp_probs = (1 - p(i))*cmp_probs;
			zicmp_probs(0) += p(i);
			x(i) = qdiscrete(exp(logq(i)), zicmp_probs);
		}
	} else {
		Rcpp::stop("lambda, nu, and p must have length n");
	}

	return x;
}
