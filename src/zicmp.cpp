#include "zicmp.h"
#include "cmp.h"
#include "z.h"
#include "util.h"

Rcpp::NumericVector d_zicmp(const Rcpp::NumericVector& x, double lambda,
	double nu, double p, bool take_log, double hybrid_tol,
	double truncate_tol, double ymax)
{
	unsigned int n = x.size();

	// Normalizing constant of CMP (not ZICMP)
	double lnormconst = z_hybrid(lambda, nu, true, hybrid_tol, truncate_tol, ymax);

	// TBD: We need to check for integer values, here and probably in
	// d_cmp. Or we can take integers as input...
	const Rcpp::LogicalVector& ind_zero = (x == 0);

	Rcpp::NumericVector out(n);
	for (unsigned int i = 0; i < n; i++) {
		// Use the property: log(a + b) = log(a) + log(1 + b/a); here, a
		// represents the count part of the density and b represents the
		// zero-inflated part.
		double log_a = log(1-p) + x(i)*log(lambda) - nu*lgamma(x(i) + 1) - lnormconst;
		double log_b = ind_zero(i) * log(p);
		double lp = logadd(log_a, log_b);
		out(i) = lp;
	}

	if (take_log) {
		return out;
	} else {
		return exp(out);
	}
}

double loglik_zicmp(const Rcpp::NumericVector& x,
	const Rcpp::NumericVector& lambda, const Rcpp::NumericVector& nu,
	const Rcpp::NumericVector& p, double hybrid_tol, double truncate_tol,
	double ymax)
{
	unsigned int n = x.size();
	double out = 0;

	for (unsigned int i = 0; i < n; i++) {
		Rcpp::NumericVector lp_vec = d_zicmp(Rcpp::NumericVector::create(x(i)),
			lambda(i), nu(i), p(i), true, hybrid_tol, truncate_tol, ymax);
		out += lp_vec(0);
	}

	return out;
}

Rcpp::NumericVector q_zicmp(const Rcpp::NumericVector& logq, double lambda,
	double nu, double p, double hybrid_tol, double truncate_tol, double ymax)
{
	// Normalizing constant of CMP (not ZICMP)
	double lnormconst = z_hybrid(lambda, nu, true, hybrid_tol, truncate_tol, ymax);

	double logq_max = Rcpp::max(logq);
	std::vector<double> all_lcp_vec;

	// Compute all of the probabilities we'll need, on the log-scale

	// The expression for f(0)
	// Use the property: log(a + b) = log(a) + log(1 + b/a).
	double log_a = log(1-p) - lnormconst;
	double log_b = log(p);
	double lp = logadd(log_a, log_b);
	double lcp = lp;
	all_lcp_vec.push_back(lcp);

	for (unsigned int j = 1; j <= ymax; j++) {
		// Do summation on the log-scale.
		// Use the property: log(a + b) = log(a) + log(1 + b/a).
		lp = j*log(lambda) - nu*lgamma(j+1) + log(1-p) - lnormconst;
		lcp = logadd(lcp, lp);
		all_lcp_vec.push_back(lcp);

		if (lcp >= logq_max) {
			break;
		}
	}

	Rcpp::NumericVector all_lcp(all_lcp_vec.begin(), all_lcp_vec.end());

	unsigned int n = logq.size();
	Rcpp::NumericVector out(n);
	for (unsigned int i = 0; i < n; i++) {
		out(i) = qdiscrete(logq(i), all_lcp);
	}

	return out;
}
