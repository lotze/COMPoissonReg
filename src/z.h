#ifndef Z_H
#define Z_H

#include <Rcpp.h>

double z_trunc(double lambda, double nu, double tol, bool take_log, double ymax);

double z_approx(double lambda, double nu, bool take_log);

double z_hybrid(double lambda, double nu, bool take_log, double hybrid_tol,
	double truncate_tol, double ymax);

// Brute force computation of the z-function. This computation is done at the
// original scale (as opposed to the log-scale), so it becomes unstable when
// the magnitudes of the terms become very large.
// [[Rcpp::export]]
Rcpp::NumericVector z_trunc(const Rcpp::NumericVector& lambda,
	const Rcpp::NumericVector& nu, double tol, bool take_log,
	double ymax);

// Approximation to the z-function from Shmueli et al (JRSS-C, 2005)
// [[Rcpp::export]]
Rcpp::NumericVector z_approx(const Rcpp::NumericVector& lambda,
	const Rcpp::NumericVector& nu, bool take_log);

// Sum (j>=0) { lambda^j / (j!)^nu }
// If lambda^{-1/nu} is small, it is more efficient (and often numerically more stable)
// to use an approximation at the log-scale.
// [[Rcpp::export]]
Rcpp::NumericVector z_hybrid(const Rcpp::NumericVector& lambda,
	const Rcpp::NumericVector& nu, bool take_log,
	double hybrid_tol, double truncate_tol, double ymax);

#endif
