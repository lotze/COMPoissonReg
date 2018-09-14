#ifndef UTIL_H
#define UTIL_H

#include <Rcpp.h>

Rcpp::IntegerVector which(const Rcpp::LogicalVector& x);

// [[Rcpp::export]]
double qdiscrete(double q, const Rcpp::NumericVector& p);

#endif