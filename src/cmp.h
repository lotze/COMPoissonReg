#ifndef CMP_H
#define CMP_H

#include <Rcpp.h>

//' Loglikelihood function for CMP.
//' 
//' @param x Vector of observed counts.
//' @param lambda Vector of lambda parameters.
//' @param nu Vector of nu parameters.
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
double loglik_cmp(const Rcpp::NumericVector& x,
	const Rcpp::NumericVector& lambda, const Rcpp::NumericVector& nu,
	double hybrid_tol, double truncate_tol, double ymax);

//' Density function for CMP.
//' 
//' @param x Vector of density arguments.
//' @param lambda Rate parameter.
//' @param nu Dispersion parameter.
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
Rcpp::NumericVector d_cmp(const Rcpp::NumericVector& x, double lambda,
	double nu, bool take_log, bool normalize, double hybrid_tol,
	double truncate_tol, double ymax);

//' CDF function for CMP.
//' 
//' @param x Vector of CDF arguments.
//' @param lambda Rate parameter.
//' @param nu Dispersion parameter.
//' @param take_log Return CDF value on the log-scale if \code{TRUE}. Otherwise,
//' return value on the original scale.
//' @param hybrid_tol Tolerance for truncation.
//' @param truncate_tol Tolerance for when to use approximation vs. truncation.
//' @param ymax Maximum value of y to consider.
//' 
//' @return CDF values.
//' 
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector p_cmp(const Rcpp::NumericVector& x, double lambda,
	double nu, double hybrid_tol, double truncate_tol, double ymax);


//' Quantile function for CMP.
//' 
//' @param logq Vector of probabilities on the log scale.
//' @param lambda Rate parameter.
//' @param nu Dispersion parameter.
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
Rcpp::NumericVector q_cmp(const Rcpp::NumericVector& logq, double lambda,
	double nu, double hybrid_tol, double truncate_tol, double ymax);

//' Draw variates from CMP.
//' 
//' @param n Number of draws to generate.
//' @param lambda Rate parameter.
//' @param nu Dispersion parameter.
//' @param hybrid_tol Tolerance for truncation.
//' @param truncate_tol Tolerance for when to use approximation vs. truncation.
//' @param ymax Maximum value of y to consider.
//' 
//' @return Draws.
//' 
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector r_cmp(unsigned int n, double lambda, double nu,
	double hybrid_tol, double truncate_tol, double ymax);

#endif
