#ifndef UTIL_H
#define UTIL_H

#include <Rcpp.h>

Rcpp::IntegerVector which(const Rcpp::LogicalVector& x);

// [[Rcpp::export]]
double qdiscrete(double q, const Rcpp::NumericVector& p);

// [[Rcpp::export]]
double qdiscrete2(double logq, const Rcpp::NumericVector& logp);

// [[Rcpp::export]]
double logsumprobs(const Rcpp::NumericVector& logprob, unsigned int baseidx);

// [[Rcpp::export]]
Rcpp::NumericVector logcumprobs(const Rcpp::NumericVector& logprob);

#endif
