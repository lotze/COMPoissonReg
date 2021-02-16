#ifndef UTIL_H
#define UTIL_H

#include <Rcpp.h>

double log(unsigned int x);

Rcpp::IntegerVector which(const Rcpp::LogicalVector& x);

// [[Rcpp::export]]
Rcpp::NumericVector linspace(double start, double end, unsigned int N);

//' @export
// [[Rcpp::export]]
unsigned int qdiscrete(double q, const Rcpp::NumericVector& cp);

// [[Rcpp::export]]
double logsumprobs(const Rcpp::NumericVector& logprob, unsigned int baseidx);

// [[Rcpp::export]]
Rcpp::NumericVector logcumprobs(const Rcpp::NumericVector& logprob);

#endif
