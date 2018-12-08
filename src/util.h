#ifndef UTIL_H
#define UTIL_H

#include <Rcpp.h>

double log(unsigned int x);

Rcpp::IntegerVector which(const Rcpp::LogicalVector& x);

// [[Rcpp::export]]
double qdiscrete(double q, const Rcpp::NumericVector& p, bool log_scale = false);

// [[Rcpp::export]]
double logsumprobs(const Rcpp::NumericVector& logprob, unsigned int baseidx);

// [[Rcpp::export]]
Rcpp::NumericVector logcumprobs(const Rcpp::NumericVector& logprob);

#endif
