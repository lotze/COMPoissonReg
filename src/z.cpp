#include "z.h"
#include "util.h"

std::pair<double, unsigned int> truncate(double lambda, double nu, double tol,
	double ymax)
{
	double log_tol = log(tol);
	double diff = R_PosInf;
	unsigned int y;

	// Compute for y = 0 outside of the loop to avoid errors in log_z_trunc sum
	double lp = -nu*lgamma(1);
	double log_z_trunc = lp;

	for (y = 1; (diff > log_tol) && (y < ymax); y++) {
		lp = y*log(lambda) - nu*lgamma(y + 1);

		// Sum normalizing constant on the log scale. Try to avoid too much
		// numerical error. Use the property: log(a + b) = log(a) + log(1 + b/a)
		log_z_trunc = logadd(log_z_trunc, lp);

		// log_ratio needs to be < 0 before we can use our Stirling
		// approximation to bound the sum. Until then, consider diff to
		// be infinity.
		double log_ratio = log(lambda) + nu - nu*log(y+1);
		if (log_ratio < 0) {
			double log_delta = -nu / 2 * log(2 * M_PI) -
				nu * (y + 3 / 2.0) * log(y + 1) +
				(y + 1) * (nu + log(lambda)) -
				log1p(-lambda * exp(nu) / pow(y + 1, nu));
			diff = log_delta - log_z_trunc;
		}

	 	if (y % 10000 == 0) {
	 		R_CheckUserInterrupt();
	 	}
	}

	if (std::isinf(diff)) {
		char msg[512];
		sprintf(msg,
			"Terms of normalizing constant CMP(%g, %g) could not be bounded by "
			"a geometric series when y <= %g. Consider adjusting the controls "
			"ymax, hybrid.tol, and truncate.tol",
			lambda, nu, ymax);
		Rf_warning(msg);
	} else if (diff > log_tol) {
		char msg[512];
		sprintf(msg,
			"Absolute relative error %g was larger than tolerance %g with "
			"CMP(%g, %g) truncated to %g. Consider adjusting the controls "
			"ymax, hybrid.tol, and COMPoissonReg.truncate.tol",
			exp(diff), tol, lambda, nu, ymax);
		Rf_warning(msg);
	}

	return std::pair<double, unsigned int>(log_z_trunc, y);
}

unsigned int y_trunc(double lambda, double nu, double tol, double ymax)
{
	const std::pair<double, unsigned int>& ret_pair = truncate(lambda, nu, tol, ymax);
	return ret_pair.second;
}

double z_trunc(double lambda, double nu, double tol, bool take_log, double ymax)
{
	const std::pair<double, unsigned int>& ret_pair = truncate(lambda, nu, tol, ymax);
	double log_z_trunc = ret_pair.first;
	if (take_log) {
		return log_z_trunc;
	} else {
		return exp(log_z_trunc);
	}
}

double z_approx(double lambda, double nu, bool take_log)
{
	double out = nu*exp(1/nu * log(lambda)) -
		(nu-1)/(2*nu) * log(lambda) -
		(nu-1)/2 * log(2*M_PI) - 0.5*log(nu);

	if (take_log) {
		return out;
	} else {
		return exp(out);
	}
}

double z_hybrid(double lambda, double nu, bool take_log, double hybrid_tol,
	double truncate_tol, double ymax)
{
	// Use the asymptotic approximation if lambda^(-1/nu) < hybrid_tol
	bool use_approx = -1/nu * log(lambda) < log(hybrid_tol);
	double out;
	
	if (use_approx) {
		out = z_approx(lambda, nu, take_log);
	} else {
		out = z_trunc(lambda, nu, truncate_tol, take_log, ymax);
	}

	return out;
}

Rcpp::NumericVector z_trunc(const Rcpp::NumericVector& lambda,
	const Rcpp::NumericVector& nu, double tol, bool take_log,
	double ymax)
{
	unsigned int n = lambda.size();
	if (n != nu.size()) {
		Rcpp::stop("Length of lambda must be same as length of nu");
	}

	Rcpp::NumericVector out(n);
	for (unsigned int i = 0; i < n; i++) {
		out(i) = z_trunc(lambda(i), nu(i), tol, take_log, ymax);
	}

	return out;
}

Rcpp::NumericVector z_approx(const Rcpp::NumericVector& lambda,
	const Rcpp::NumericVector& nu, bool take_log)
{
	unsigned int n = lambda.size();
	if (n != nu.size()) {
		Rcpp::stop("Length of lambda must be same as length of nu");
	}

	Rcpp::NumericVector out(n);
	for (unsigned int i = 0; i < n; i++) {
		out(i) = z_approx(lambda(i), nu(i), take_log);
	}

	return out;
}

Rcpp::NumericVector z_hybrid(const Rcpp::NumericVector& lambda,
	const Rcpp::NumericVector& nu, bool take_log,
	double hybrid_tol, double truncate_tol, double ymax)
{
	unsigned int n = lambda.size();
	if (n != nu.size()) {
		Rcpp::stop("lambda and nu must be the same length");
	}

	Rcpp::NumericVector out(n);
	for (unsigned int i = 0; i < n; i++) {
		out(i) = z_hybrid(lambda(i), nu(i), take_log, hybrid_tol,
			truncate_tol, ymax);
	}

	return out;
}

Rcpp::IntegerVector y_trunc(const Rcpp::NumericVector& lambda,
	const Rcpp::NumericVector& nu, double tol, double ymax)
{
	unsigned int n = lambda.size();
	if (n != nu.size()) {
		Rcpp::stop("Length of lambda must be same as length of nu");
	}

	Rcpp::IntegerVector out(n);
	for (unsigned int i = 0; i < n; i++) {
		out(i) = y_trunc(lambda(i), nu(i), tol, ymax);
	}

	return out;
}

