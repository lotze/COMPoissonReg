#ifndef CMP_H
#define CMP_H

#include <Rcpp.h>

// TBD: I think we might want to delete this, but not sure if we can yet.
// @export
// // [[Rcpp::export]]
//Rcpp::NumericVector allprobs_cmp(double lambda, double nu, double tol,
//	bool take_log, double ymax, bool normalize);

// [[Rcpp::export]]
double loglik_cmp(const Rcpp::NumericVector& x,
	const Rcpp::NumericVector& lambda, const Rcpp::NumericVector& nu,
	double hybrid_tol, double truncate_tol, double ymax);

// [[Rcpp::export]]
Rcpp::NumericVector d_cmp(const Rcpp::NumericVector& x, double lambda,
	double nu, bool take_log, bool normalize, double hybrid_tol,
	double truncate_tol, double ymax);

// [[Rcpp::export]]
Rcpp::NumericVector p_cmp(const Rcpp::NumericVector& x, double lambda,
	double nu, double hybrid_tol, double truncate_tol, double ymax);

// Assume that quantiles q are given on the log scale.
// Work on the log-scale for stability.
// [[Rcpp::export]]
Rcpp::NumericVector q_cmp(const Rcpp::NumericVector& logq, double lambda,
	double nu, double hybrid_tol, double truncate_tol, double ymax);

// Produce n iid draws.
// If n is a scalar, take n(i) to be n
// [[Rcpp::export]]
Rcpp::NumericVector r_cmp(unsigned int n, double lambda, double nu,
	double hybrid_tol, double truncate_tol, double ymax);

#endif
