#include <Rcpp.h>
#include <list>
#include "util.h"

// Solaris gives errors on CRAN if we do not define this.
double log(unsigned int x)
{
	return log(double(x));
}

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
// that p(k) = 1 - p(0) - ... - p(k-1). If log_scale = TRUE, interpret q and p
// as log-probabilities
double qdiscrete(double q, const Rcpp::NumericVector& p, bool log_scale)
{
	unsigned int k = p.size();
	Rcpp::IntegerVector idx;

	if (log_scale) {
		Rcpp::NumericVector lcp = logcumprobs(p);
		idx = which(q < lcp);
	} else {
		Rcpp::NumericVector cumprobs = Rcpp::cumsum(p);
		idx = which(q < cumprobs);
	}

	if (idx.size() > 0) {
		return idx(0);
	} else {
		return k;
	}
}

// Compute log(p[0] + ... + p[k-1]) given log(p[0]), ..., log(p[k-1]).
// Treat p[baseidx] as the "normalizing" probability. Uses the method from
// https://en.wikipedia.org/wiki/List_of_logarithmic_identities#Summation/subtraction
// for stable calculation on the log scale.
double logsumprobs(const Rcpp::NumericVector& logprob, unsigned int baseidx)
{
        unsigned int k = logprob.size();
		Rcpp::NumericVector logprob_rest(k-1);

		for (unsigned int i = 0, j = 0; i < k; i++) {
			if (i != baseidx) {
				logprob_rest(j) = logprob(i);
				j++;
			}
		}

        return logprob(baseidx) + log1p(Rcpp::sum(Rcpp::exp(logprob_rest - logprob(baseidx))));
}

// Compute cumulative probabilities on the log scale. Uses the method from
// https://en.wikipedia.org/wiki/List_of_logarithmic_identities#Summation/subtraction
// for stable calculation on the log scale.
Rcpp::NumericVector logcumprobs(const Rcpp::NumericVector& logprob)
{
        unsigned int k = logprob.size();
		Rcpp::NumericVector logcp(k);

		logcp(0) = logprob(0);
		for (unsigned int i = 1; i < k; i++) {
			logcp(i) = logcp(i-1) + log1p(exp(logprob(i) - logcp(i-1)));
		}

        return logcp;
}
