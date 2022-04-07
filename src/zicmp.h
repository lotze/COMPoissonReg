#ifndef ZICMP_H
#define ZICMP_H

#include <Rcpp.h>

/*
* Compute ZICMP density. See manual entry for the exported version of this
* function below.
*/
double d_zicmp(unsigned int x, double lambda, double nu, double p,
	bool take_log, double hybrid_tol, double truncate_tol, double ymax);

//' Loglikelihood function for ZICMP.
//' 
//' @param x Vector of observed counts.
//' @param lambda Vector of lambda parameters.
//' @param nu Vector of nu parameters.
//' @param p Vector of p parameters.
//' @param hybrid_tol Tolerance for truncation.
//' @param truncate_tol Tolerance for when to use approximation vs. truncation.
//' @param ymax Maximum value of y to consider.
//' 
//' @details
//' The vectors \code{x}, \code{lambda}, and \code{nu} must ne the same length.
//' 
//' @return Value of loglikelihood.
//' 
//' @noRd
// [[Rcpp::export]]
double loglik_zicmp(const Rcpp::IntegerVector& x,
	const Rcpp::NumericVector& lambda, const Rcpp::NumericVector& nu,
	const Rcpp::NumericVector& p, double hybrid_tol, double truncate_tol,
	double ymax);

//' Density function for ZICMP.
//' 
//' @param x Vector of density arguments.
//' @param lambda Rate parameter.
//' @param nu Dispersion parameter.
//' @param p Zero-inflation parameter.
//' @param take_log Return density on the log-scale if \code{TRUE}. Otherwise,
//' return density on the original scale.
//' @param hybrid_tol Tolerance for truncation.
//' @param truncate_tol Tolerance for when to use approximation vs. truncation.
//' @param ymax Maximum value of y to consider.
//' 
//' @return Density values.
//' 
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector d_zicmp(const Rcpp::IntegerVector& x, double lambda,
	double nu, double p, bool take_log, double hybrid_tol,
	double truncate_tol, double ymax);

//' Quantile function for ZICMP.
//' 
//' @param logq Vector of probabilities on the log scale.
//' @param lambda Rate parameter.
//' @param nu Dispersion parameter.
//' @param p Zero-inflation parameter.
//' @param hybrid_tol Tolerance for truncation.
//' @param truncate_tol Tolerance for when to use approximation vs. truncation.
//' @param ymax Maximum value of y to consider.
//' 
//' @details
//' Operates on the log-scale for stability.
//' 
//' @return Quantile values.
//' 
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector q_zicmp(const Rcpp::NumericVector& logq,
	double lambda, double nu, double p, double hybrid_tol, double truncate_tol,
	double ymax);

#endif
