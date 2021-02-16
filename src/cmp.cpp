#include "cmp.h"
#include <algorithm>
#include <vector>
#include "z.h"
#include "util.h"

/*
// TBD: Do we still need this, or can it be removed?
// Enumerate the terms lambda^y / (y!)^nu for y >= 0, until they become small,
// or until y = ymax is reached.
Rcpp::NumericVector allprobs_cmp(double lambda, double nu, double tol,
	bool take_log, double ymax, bool normalize)
{
	std::vector<double> logp_unnorm;
	double log_tol = log(tol);
	double diff = R_PosInf;
	unsigned int y;

	// Compute for y = 0 outside of the loop to avoid errors in log_Z_trunc sum
	double lp = -nu*lgamma(1);
	double log_Z_trunc = lp;
	logp_unnorm.push_back(lp);

	for (y = 1; (diff > log_tol) && (y < ymax); y++) {
		lp = y*log(lambda) - nu*lgamma(y + 1);
		logp_unnorm.push_back(lp);

		// Sum normalizing constant on the log scale. Try to avoid too much numerical error.
		// Using the property: log(a + b) = log(a) + log(1 + b/a)
		log_Z_trunc += log1p(exp(lp - log_Z_trunc));

		double log_ratio = log(lambda) + nu - nu*log(y+1);
		if (log_ratio < 0) {
			double log_Delta = -nu/2 * log(2*M_PI) - nu*(y + 3/2.0) * log(y + 1) +
				(y + 1)*(nu + log(lambda)) - log1p(-lambda * exp(nu) / pow(y + 1, nu));
			diff = log_Delta - log_Z_trunc;
		}

	 	if (y % 10000 == 0) {
	 		R_CheckUserInterrupt();
	 	}
	}

	if (y == ymax) {
		char msg[128];
		sprintf(msg, "%s\n\toptions(COMPoissonReg.ymax = %g)\n",
			"Larger numbers may be needed for CMP. Try increasing this setting:", ymax);
		Rf_warning(msg);
	}

	Rcpp::NumericVector logp(logp_unnorm.begin(), logp_unnorm.end());

	if (normalize) {
		logp = logp - log_Z_trunc;
	}

	if (take_log) {
		return logp;
	} else {
		return exp(logp);
	}
}
*/

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
		out = out - z_hybrid(lambda, nu, true, hybrid_tol, truncate_tol, ymax);
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

	double lnormconst = z_hybrid(lambda, nu, true, hybrid_tol,
		truncate_tol, ymax);

	unsigned int x_max = int(std::min(double(Rcpp::max(x)), ymax));
	const Rcpp::NumericVector& x_range = linspace(0, x_max, x_max + 1);
	const Rcpp::NumericVector& all_lp_unnorm = d_cmp(x_range, lambda, nu, true,
			false, hybrid_tol, truncate_tol, ymax);

	for (unsigned int i = 0; i < n; i++) {
		double lp = lnormconst;
		for (unsigned int j = 0; j <= x(i) && j <= x_max; j++) {
			// Do summation on the log-scale.
			// Use the property: log(a + b) = log(a) + log(1 + b/a).
			lp += log1p(exp(all_lp_unnorm(j) - lp));
		}
		out(i) = lp;
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
		lcp += log1p(exp(lp - lcp));
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
		out(i) = qdiscrete(logq(i), all_lcp);
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
