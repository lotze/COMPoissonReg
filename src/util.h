#ifndef UTIL_H
#define UTIL_H

#include <Rcpp.h>

/* The following macro helps us to use the property:
*  log(x + y) = log(x) + log(1 + y/x)
*             = log(x) + log1p(exp(log(y) - log(x)))
*  When x and y are of very different magnitudes, this is more stable
*  when the x is taken to be the larger of the two.
*  https://en.wikipedia.org/wiki/List_of_logarithmic_identities#Summation
*/
#define logadd(logx, logy) \
	((logx + log1p(exp(logy - logx))) * (logx > logy) + \
		(logy + log1p(exp(logx - logy))) * (logy >= logx))

double log(unsigned int x);

Rcpp::IntegerVector which(const Rcpp::LogicalVector& x);

// [[Rcpp::export]]
Rcpp::NumericVector linspace(double start, double end, unsigned int N);

//' @export
// [[Rcpp::export]]
unsigned int qdiscrete(double q, const Rcpp::NumericVector& cp);

#endif
