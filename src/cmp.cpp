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
	double psi = -0.5772157;	// Euler–Mascheroni constant
	double y = 0;
	double deriv = R_PosInf;
	while (deriv > 0 && y <= ymax) {
		delta = y*log(lambda) - nu*lgamma(y+1);
		logp_unnorm.push_back(delta);
		psi += 1 / (y+1);
		deriv = log(lambda) - nu*psi;
		y++;
	}

	double log_tol = log(tol);
	while (delta > log_tol && y <= ymax) {
		delta = y*log(lambda) - nu*lgamma(y+1);
		logp_unnorm.push_back(delta);
		y++;
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
	if (n != lambda.size() || n != nu.size()) {
		Rcpp::stop("lambda and nu must have length n");
	}

	Rcpp::NumericVector fx(n);
	fx.fill(0);

	for (unsigned int i = 0; i < n; i++) {
		Rcpp::NumericVector allprobs = cmp_allprobs(lambda(i), nu(i), tol);
		if (x(i) <= allprobs.size() - 1) {
			fx(i) = allprobs(x(i));
		}
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
	if (n != lambda.size() || n != nu.size()) {
		Rcpp::stop("lambda and nu must have length n");
	}

	Rcpp::NumericVector Fx(n);
	Fx.fill(1);

	for (unsigned int i = 0; i < n; i++) {
		Rcpp::NumericVector allprobs = cmp_allprobs(lambda(i), nu(i), tol);
		Rcpp::NumericVector cumprobs = Rcpp::cumsum(allprobs);
		if (x(i) <= allprobs.size() - 1) {
			Fx(i) = cumprobs(x(i));
		}
	}

	return Fx;
}

// Assume that quantiles q are not given on the log scale
Rcpp::NumericVector qcmp_cpp(const Rcpp::NumericVector& q,
	const Rcpp::NumericVector& lambda, const Rcpp::NumericVector& nu,
	double tol)
{
	unsigned int n = q.size();
	if (n != lambda.size() || n != nu.size()) {
		Rcpp::stop("lambda and nu must have length n");
	}

	Rcpp::NumericVector x(n);
	for (unsigned int i = 0; i < n; i++) {
		// Find the q(i) quantile on the log-scale
		Rcpp::NumericVector all_logprobs = cmp_allprobs(lambda(i), nu(i), tol, true);
		x(i) = qdiscrete(log(q(i)), all_logprobs, true);
	}

	return x;
}

// Produce n(i) iid draws for each lambda(i) and nu(i).
// If n is a scalar, take n(i) to be n
Rcpp::NumericVector rcmp_cpp(unsigned int n, const Rcpp::NumericVector& lambda,
	const Rcpp::NumericVector& nu, double tol)
{
	if (n != lambda.size() || n != nu.size()) {
		Rcpp::stop("lambda and nu must have length n");
	}

	Rcpp::NumericVector u = Rcpp::runif(n, 0.0, 1.0);
	return qcmp_cpp(u, lambda, nu, tol);
}
