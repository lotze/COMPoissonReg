#ifndef Z_H
#define Z_H

#include <Rcpp.h>

/*
* Compute the normalizing constant using an appropriate truncation Also return
* the truncation value instead of the normalizing constant. Other functions can
* use this to find a suitable value for truncated operations.
*/
std::pair<double, unsigned int> truncate(double lambda, double nu, double tol, double ymax);

/*
* A simple accessor to the second element of truncate's return value
*/
unsigned int y_trunc(double lambda, double nu, double tol, double ymax);

/*
* A simple accessor to the first element of truncate's return value
*/
double z_trunc(double lambda, double nu, double tol, bool take_log, double ymax);

/*
* Compute normalizing constant by asymptotic approximation
*/
double z_approx(double lambda, double nu, bool take_log);

/*
* Hybrid method to compute normalizing constant. See manual entry for the
* exported version of this function below.
*/
double z_hybrid(double lambda, double nu, bool take_log, double hybrid_tol,
	double truncate_tol, double ymax);

// TBD: Update this
// Brute force computation of the z-function. This computation is done at the
// original scale (as opposed to the log-scale), so it becomes unstable when
// the magnitudes of the terms become very large.
// [[Rcpp::export]]
Rcpp::NumericVector z_trunc(const Rcpp::NumericVector& lambda,
	const Rcpp::NumericVector& nu, double tol, bool take_log,
	double ymax);

// TBD: Update this
// Approximation to the z-function from Shmueli et al (JRSS-C, 2005)
// [[Rcpp::export]]
Rcpp::NumericVector z_approx(const Rcpp::NumericVector& lambda,
	const Rcpp::NumericVector& nu, bool take_log);

// TBD: Update this
// Sum (j>=0) { lambda^j / (j!)^nu }
// If lambda^{-1/nu} is small, it is more efficient (and often numerically more stable)
// to use an approximation at the log-scale.
// [[Rcpp::export]]
Rcpp::NumericVector z_hybrid(const Rcpp::NumericVector& lambda,
	const Rcpp::NumericVector& nu, bool take_log,
	double hybrid_tol, double truncate_tol, double ymax);

// TBD: Update this
// Return the truncation value computed in z_trunc. This can be used by other
// functions to compute truncated expressions.
// [[Rcpp::export]]
Rcpp::IntegerVector y_trunc(const Rcpp::NumericVector& lambda,
	const Rcpp::NumericVector& nu, double tol, double ymax);

#endif
