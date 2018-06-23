#include <Rcpp.h>
#include <list>
#include "COMPoissonReg.h"

// A version of the "which" function, since it does not appear to be
// provided in Rcpp Sugar yet.
Rcpp::IntegerVector which(const Rcpp::LogicalVector& x)
{
	std::list<int> idx;
	for (int i = 0; i < x.size(); i++) {
		if (x(i)) {
			idx.push_back(i);
		}
	}

	return Rcpp::IntegerVector(idx.begin(), idx.end());
}

// Quantile function for a discrete distribution with values (0, 1, ..., k-1)
// and probabilities p = (p(0), ..., p(k-1)). If p does not sum to 1, we assume
// that p(k) = 1 - p(0) - ... - p(k-1).
double qdiscrete(double q, const Rcpp::NumericVector& p)
{
	unsigned int k = p.size();
	Rcpp::NumericVector cumprobs = Rcpp::cumsum(p);

	Rcpp::IntegerVector idx = which(q < cumprobs);
	if (idx.size() > 0) {
		return idx(0);
	} else {
		return k;
	}
}
