#include "zicmp.h"
#include "cmp.h"
#include "z.h"
#include "util.h"

double d_zicmp(unsigned int x, double lambda, double nu, double p,
	bool take_log, double hybrid_tol, double truncate_tol, double ymax)
{
	// Normalizing constant of CMP (not ZICMP)
	// Here it makes less sense to return unnormalized density values, so
	// we don't provide it as an option.
	double lnormconst = z_hybrid(lambda, nu, true, hybrid_tol, truncate_tol, ymax);

	// Check for zero
	bool ind_zero = (x == 0);

	// Use the property: log(a + b) = log(a) + log(1 + b/a); here, a
	// represents the count part of the density and b represents the
	// zero-inflated part.
	double log_a = log(1-p) + x*log(lambda) - nu*lgamma(x+1) - lnormconst;
	double log_b = log(ind_zero * p);
	double out = logadd(log_a, log_b);

	if (take_log) {
		return out;
	} else {
		return exp(out);
	}
}


Rcpp::NumericVector d_zicmp(const Rcpp::IntegerVector& x, double lambda,
	double nu, double p, bool take_log, double hybrid_tol,
	double truncate_tol, double ymax)
{
	unsigned int n = x.size();

	// Normalizing constant of CMP (not ZICMP)
	// Avoid recomputing this for each element of x
	double lnormconst = z_hybrid(lambda, nu, true, hybrid_tol, truncate_tol, ymax);

	// Check for zeros
	const Rcpp::LogicalVector& ind_zero = (x == 0);

	Rcpp::NumericVector out(n);
	for (unsigned int i = 0; i < n; i++) {
		// Use the property: log(a + b) = log(a) + log(1 + b/a); here, a
		// represents the count part of the density and b represents the
		// zero-inflated part.
		double log_a = log(1-p) + x(i)*log(lambda) - nu*lgamma(x(i) + 1) - lnormconst;
		double log_b = log(ind_zero(i) * p);
		double lp = logadd(log_a, log_b);
		out(i) = lp;
	}

	if (take_log) {
		return out;
	} else {
		return exp(out);
	}
}

double loglik_zicmp(const Rcpp::IntegerVector& x,
	const Rcpp::NumericVector& lambda, const Rcpp::NumericVector& nu,
	const Rcpp::NumericVector& p, double hybrid_tol, double truncate_tol,
	double ymax)
{
	unsigned int n = x.size();
	double out = 0;

	for (unsigned int i = 0; i < n; i++) {
		out += d_zicmp(x(i), lambda(i), nu(i), p(i), true, hybrid_tol, truncate_tol, ymax);
	}

	return out;
}

Rcpp::NumericVector q_zicmp(const Rcpp::NumericVector& logq, double lambda,
	double nu, double p, double hybrid_tol, double truncate_tol, double ymax)
{
	// Normalizing constant of CMP (not ZICMP)
	// Since we're using the truncated method below, we'll compute the normcost by
	// truncation too. We can use the same call to get our upper truncation bound.
	const std::pair<double, unsigned int>& ret_pair = truncate(lambda, nu,
		truncate_tol, ymax);
	double lnormconst = ret_pair.first;
	unsigned int M = ret_pair.second;

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

	for (unsigned int j = 1; j <= M; j++) {
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
		out(i) = q_discrete(logq(i), all_lcp);
	}

	return out;
}
