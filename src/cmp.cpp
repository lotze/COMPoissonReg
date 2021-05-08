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

	double lnormconst = z_hybrid(lambda, nu, true, hybrid_tol, truncate_tol, ymax);

	unsigned int x_max = int(std::min(double(Rcpp::max(x)), ymax));

	for (unsigned int i = 0; i < n; i++) {
		double lcp = -lnormconst;
		for (unsigned int j = 1; j <= x(i) && j <= x_max; j++) {
			// Do summation on the log-scale.
			// Use the property: log(a + b) = log(a) + log(1 + b/a).
			double lp = j*log(lambda) - nu*lgamma(j + 1) - lnormconst;
			lcp = logadd(lcp, lp);
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
	double lnormconst = z_hybrid(lambda, nu, true, hybrid_tol,
		truncate_tol, ymax);
	
	double logq_max = Rcpp::max(logq);
	std::vector<double> all_lcp_vec;

	// Compute all of the probabilities we'll need, on the log-scale
	// Initialize with f(0)
	double lcp = -lnormconst;
	all_lcp_vec.push_back(lcp);

	for (unsigned int j = 1; j <= ymax; j++) {
		// Do summation on the log-scale.
		// Use the property: log(a + b) = log(a) + log(1 + b/a).
		double lp = j*log(lambda) - nu*lgamma(j+1) - lnormconst;
		lcp = logadd(lcp, lp);
		all_lcp_vec.push_back(lcp);

	 	if (j % 10000 == 0 && j > 0) {
	 		R_CheckUserInterrupt();
	 	}

	 	if (lcp > logq_max) {
	 		break;
	 	}
	}

	Rcpp::NumericVector all_lcp(all_lcp_vec.begin(), all_lcp_vec.end());

	// Rcpp::print(all_lcp);

	unsigned int n = logq.size();
	Rcpp::NumericVector out(n);
	for (unsigned int i = 0; i < n; i++) {
		out(i) = q_discrete(logq(i), all_lcp);
		// Rprintf("Checkpoint 1: We drew %g for q(%d) = %g\n", out(i), i, exp(logq(i)));
		R_CheckUserInterrupt();
	}

	return out;
}

Rcpp::NumericVector r_cmp(unsigned int n, double lambda, double nu,
	double hybrid_tol, double truncate_tol, double ymax)
{
	const Rcpp::NumericVector& u = Rcpp::runif(n, 0.0, 1.0);
	return q_cmp(Rcpp::log(u), lambda, nu, hybrid_tol, truncate_tol, ymax);
}
