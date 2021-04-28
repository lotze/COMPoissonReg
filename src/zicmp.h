#ifndef ZICMP_H
#define ZICMP_H

#include <Rcpp.h>

// [[Rcpp::export]]
double loglik_zicmp(const Rcpp::NumericVector& x,
	const Rcpp::NumericVector& lambda, const Rcpp::NumericVector& nu,
	const Rcpp::NumericVector& p, double hybrid_tol, double truncate_tol,
	double ymax);

// [[Rcpp::export]]
Rcpp::NumericVector d_zicmp(const Rcpp::NumericVector& x, double lambda,
	double nu, double p, bool take_log, double hybrid_tol,
	double truncate_tol, double ymax);

// Assume that quantiles q are given on the log scale
// [[Rcpp::export]]
Rcpp::NumericVector q_zicmp(const Rcpp::NumericVector& logq,
	double lambda, double nu, double p, double hybrid_tol, double truncate_tol,
	double ymax);

#endif
