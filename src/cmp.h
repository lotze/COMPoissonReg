#ifndef CMP_H
#define CMP_H

#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::NumericVector cmp_allprobs(double lambda, double nu, double tol = 1e-6,
	bool take_log = false, double ymax = 1e100, bool normalize = true);

// [[Rcpp::export]]
Rcpp::NumericVector dcmp_cpp(const Rcpp::NumericVector& x,
	const Rcpp::NumericVector& lambda, const Rcpp::NumericVector& nu,
	double tol = 1e-6, bool take_log = false, double ymax = 1e100);

// [[Rcpp::export]]
Rcpp::NumericVector pcmp_cpp(const Rcpp::NumericVector& x,
	const Rcpp::NumericVector& lambda, const Rcpp::NumericVector& nu,
	double tol = 1e-6, double ymax = 1e100);

// [[Rcpp::export]]
Rcpp::NumericVector qcmp_cpp(const Rcpp::NumericVector& logq,
	const Rcpp::NumericVector& lambda, const Rcpp::NumericVector& nu,
	double tol = 1e-6, double ymax = 1e100);

// [[Rcpp::export]]
Rcpp::NumericVector rcmp_cpp(unsigned int n, const Rcpp::NumericVector& lambda,
	const Rcpp::NumericVector& nu, double tol = 1e-6, double ymax = 1e100);

#endif
