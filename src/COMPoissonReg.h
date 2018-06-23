#ifndef COMPOISSONREG_H
#define COMPOISSONREG_H

#include <Rcpp.h>

Rcpp::IntegerVector which(const Rcpp::LogicalVector& x);

// [[Rcpp::export]]
double qdiscrete(double q, const Rcpp::NumericVector& p);

// [[Rcpp::export]]
Rcpp::NumericVector cmp_allprobs(double lambda, double nu, double tol = 1e-6,
	bool take_log = false, double ymax = 1e100, bool normalize = true);

#endif