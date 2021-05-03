#ifndef CMP_H
#define CMP_H

#include <Rcpp.h>

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
// [[Rcpp::export]]
Rcpp::NumericVector r_cmp(unsigned int n, double lambda, double nu,
	double hybrid_tol, double truncate_tol, double ymax);

#endif
