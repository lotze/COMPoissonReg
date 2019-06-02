#include <Rcpp.h>
#include <list>
#include "cmp.h"
#include "util.h"

// Enumerate the terms lambda^y / (y!)^nu for y >= 0, until they become small,
// or until y = ymax is reached.
Rcpp::NumericVector cmp_allprobs(double lambda, double nu, double tol,
	bool take_log, double ymax, bool normalize)
{
	std::list<double> logp_unnorm;

	double delta;
	double psi = -0.577215664901532;	// Euler–Mascheroni constant
	double y = 0;
	double deriv = R_PosInf;
	while (deriv > 0 && y <= ymax) {
		delta = y*log(lambda) - nu*lgamma(y+1);
		logp_unnorm.push_back(delta);
		psi += 1 / (y+1);
		deriv = log(lambda) - nu*psi;
		y++;

		if (int(y+1) % 10000 == 0) {
			R_CheckUserInterrupt();
		}
	}

	double log_tol = log(tol);
	while (delta > log_tol && y <= ymax) {
		delta = y*log(lambda) - nu*lgamma(y+1);
		logp_unnorm.push_back(delta);
		y++;

		if (int(y+1) % 10000 == 0) {
			R_CheckUserInterrupt();
		}
	}

	Rcpp::NumericVector logp(logp_unnorm.begin(), logp_unnorm.end());

	// A somewhat stable way to normalize log-probabilities
	if (normalize) {
		unsigned int idx_max = Rcpp::which_max(logp);
		logp = logp - logsumprobs(logp, idx_max);
	}

	if (take_log) {
		return logp;
	} else {
		return exp(logp);
	}
}

Rcpp::NumericVector dcmp_cpp(const Rcpp::NumericVector& x,
	const Rcpp::NumericVector& lambda, const Rcpp::NumericVector& nu,
	double tol, bool take_log)
{
	unsigned int n = x.size();
	Rcpp::NumericVector fx(n);
	fx.fill(0);

	if (lambda.size() == 1 && nu.size() == 1) {
		Rcpp::NumericVector allprobs = cmp_allprobs(lambda(0), nu(0), tol);
		for (unsigned int i = 0; i < n; i++) {
			if (x(i) <= allprobs.size() - 1) {
				fx(i) = allprobs(x(i));
			}
		}		
	} else if (lambda.size() == n || nu.size() == n) {
		for (unsigned int i = 0; i < n; i++) {
			Rcpp::NumericVector allprobs = cmp_allprobs(lambda(i), nu(i), tol);
			if (x(i) <= allprobs.size() - 1) {
				fx(i) = allprobs(x(i));
			}
		}
	} else {
		Rcpp::stop("lambda and nu must both have length n or 1");
	}

	if (take_log) {
		return log(fx);
	} else {
		return fx;
	}
}

Rcpp::NumericVector pcmp_cpp(const Rcpp::NumericVector& x,
	const Rcpp::NumericVector& lambda, const Rcpp::NumericVector& nu,
	double tol)
{
	unsigned int n = x.size();
	Rcpp::NumericVector Fx(n);
	Fx.fill(0);
	
	if (lambda.size() == 1 && nu.size() == 1) {
		Rcpp::NumericVector all_logprobs = cmp_allprobs(lambda(0), nu(0), tol, true);
		Rcpp::NumericVector lcp = logcumprobs(all_logprobs);
		for (unsigned int i = 0; i < n; i++) {
			if (x(i) <= all_logprobs.size() - 1) {
				Fx(i) = lcp(x(i));
			}
		}
	} else if (lambda.size() == n && nu.size() == n) {
		for (unsigned int i = 0; i < n; i++) {
			Rcpp::NumericVector all_logprobs = cmp_allprobs(lambda(i), nu(i), tol, true);
			Rcpp::NumericVector lcp = logcumprobs(all_logprobs);
			if (x(i) <= all_logprobs.size() - 1) {
				Fx(i) = lcp(x(i));
			}
		}
	} else {
		Rcpp::stop("lambda and nu must both have length n or 1");
	}

	return exp(Fx);
}

// Assume that quantiles q are given on the log scale
// Work on the log-scale for stability
Rcpp::NumericVector qcmp_cpp(const Rcpp::NumericVector& logq,
	const Rcpp::NumericVector& lambda, const Rcpp::NumericVector& nu,
	double tol)
{
	unsigned int n = logq.size();
	Rcpp::NumericVector x(n);
	
	if (lambda.size() == 1 && nu.size() == 1) {
		Rcpp::NumericVector all_logprobs = cmp_allprobs(lambda(0), nu(0), tol, true);	
		for (unsigned int i = 0; i < n; i++) {
			x(i) = qdiscrete(logq(i), all_logprobs, true);
		}
	} else if (lambda.size() == n && nu.size() == n) {
		for (unsigned int i = 0; i < n; i++) {
			Rcpp::NumericVector all_logprobs = cmp_allprobs(lambda(i), nu(i), tol, true);
			x(i) = qdiscrete(logq(i), all_logprobs, true);
		}
	} else {
		Rcpp::stop("lambda and nu must both have length n or 1");
	}

	return x;
}

// Produce n(i) iid draws for each lambda(i) and nu(i).
// If n is a scalar, take n(i) to be n
Rcpp::NumericVector rcmp_cpp(unsigned int n, const Rcpp::NumericVector& lambda,
	const Rcpp::NumericVector& nu, double tol)
{
	Rcpp::NumericVector u = Rcpp::runif(n, 0.0, 1.0);
	return qcmp_cpp(log(u), lambda, nu, tol);
}
