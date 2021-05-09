#include "cmp.h"
#include <algorithm>
#include <vector>
#include "z.h"
#include "util.h"

double loglik_cmp(const Rcpp::NumericVector& x,
	const Rcpp::NumericVector& lambda, const Rcpp::NumericVector& nu,
	double hybrid_tol, double truncate_tol, double ymax)
{
	unsigned int n = x.size();
	double out = 0;

	for (unsigned int i = 0; i < n; i++) {
		Rcpp::NumericVector lp_vec = d_cmp(Rcpp::NumericVector::create(x(i)),
			lambda(i), nu(i), true, true, hybrid_tol, truncate_tol, ymax);
		out += lp_vec(0);
	}

	return out;
}

Rcpp::NumericVector d_cmp(const Rcpp::NumericVector& x, double lambda, double nu,
	bool take_log, bool normalize, double hybrid_tol, double truncate_tol, double ymax)
{
	unsigned int n = x.size();

	Rcpp::NumericVector out(n);
	for (unsigned int i = 0; i < n; i++) {
		out(i) = x(i)*log(lambda) - nu*lgamma(x(i) + 1);
	}

	if (normalize) {
		double lnormcost = z_hybrid(lambda, nu, true, hybrid_tol, truncate_tol, ymax);
		out = out - lnormcost;
	}

	if (take_log) {
		return out;
	} else {
		return exp(out);
	}
}

Rcpp::NumericVector p_cmp(const Rcpp::NumericVector& x, double lambda, double nu,
	double hybrid_tol, double truncate_tol, double ymax)
{
	unsigned int n = x.size();
	Rcpp::NumericVector out(n);

	// Since we're using the truncated method below, we'll compute the normcost by
	// truncation too. We can use the same call to get our upper truncation bound.
	const std::pair<double, unsigned int>& ret_pair = truncate(lambda, nu,
		truncate_tol, ymax);
	double lnormconst = ret_pair.first;
	unsigned int M = ret_pair.second;

	unsigned int x_max = int(std::min(double(Rcpp::max(x)), double(M)));

	for (unsigned int i = 0; i < n; i++) {
		double lcp = -lnormconst;
		for (unsigned int j = 1; j <= x(i) && j <= x_max; j++) {
			// Do summation on the log-scale.
			// Use the property: log(a + b) = log(a) + log(1 + b/a).
			double lp = j*log(lambda) - nu*lgamma(j + 1) - lnormconst;
			lcp = logadd(lcp, lp);

		 	if (j % 10000 == 0) {
		 		R_CheckUserInterrupt();
	 		}
		}
		out(i) = lcp;
	}

	return Rcpp::exp(out);
}

// This method will not work for large lambda and small nu... the magnitudes
// of the numbers will become extremely large, and we won't be able to enumerate
// them.
Rcpp::NumericVector q_cmp(const Rcpp::NumericVector& logq, double lambda,
	double nu, double hybrid_tol, double truncate_tol, double ymax)
{
	
	// Since we're using the truncated method below, we'll compute the normcost by
	// truncation too. We can use the same call to get our upper truncation bound.
	const std::pair<double, unsigned int>& ret_pair = truncate(lambda, nu,
		truncate_tol, ymax);
	double lnormconst = ret_pair.first;
	unsigned int M = ret_pair.second;

	double logq_max = Rcpp::max(logq);
	std::vector<double> all_lcp_vec;

	// Compute all of the probabilities we'll need, on the log-scale.
	// Initialize with f(0).
	double lcp = -lnormconst;
	all_lcp_vec.push_back(lcp);

	for (unsigned int j = 1; j <= M; j++) {
		// Do summation on the log-scale.
		// Use the property: log(a + b) = log(a) + log(1 + b/a).
		double lp = j*log(lambda) - nu*lgamma(j+1) - lnormconst;
		lcp = logadd(lcp, lp);
		all_lcp_vec.push_back(lcp);

	 	if (j % 10000 == 0) {
	 		R_CheckUserInterrupt();
	 	}

	 	if (lcp > logq_max) {
	 		break;
	 	}
	}

	Rcpp::NumericVector all_lcp(all_lcp_vec.begin(), all_lcp_vec.end());

	unsigned int n = logq.size();
	Rcpp::NumericVector out(n);
	for (unsigned int i = 0; i < n; i++) {
		out(i) = q_discrete(logq(i), all_lcp);
	}

	return out;
}

Rcpp::NumericVector r_cmp(unsigned int n, double lambda, double nu,
	double hybrid_tol, double truncate_tol, double ymax)
{
	const Rcpp::NumericVector& u = Rcpp::runif(n, 0.0, 1.0);
	return q_cmp(Rcpp::log(u), lambda, nu, hybrid_tol, truncate_tol, ymax);
}
