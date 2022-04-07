#ifndef Z_H
#define Z_H

#include <Rcpp.h>

/*
* Compute the normalizing constant using an appropriate truncation. Also return
* the truncation value instead of the normalizing constant; other functions can
* use this to find a suitable value for truncated operations.
*/
std::pair<double, unsigned int> truncate(double lambda, double nu, double tol,
	double ymax);

/*
* A simple accessor to the second element of truncate's return value. See
* manual entry for the exported version of this function below.
*/
unsigned int y_trunc(double lambda, double nu, double tol, double ymax);

/*
* A simple accessor to the first element of truncate's return value. See manual
* entry for the exported version of this function below.
*/
double z_trunc(double lambda, double nu, double tol, bool take_log, double ymax);

/*
* Compute normalizing constant by asymptotic approximation. See manual entry
* for the exported version of this function below.
*/
double z_approx(double lambda, double nu, bool take_log);

/*
* Hybrid method to compute normalizing constant. See manual entry for the
* exported version of this function below.
*/
double z_hybrid(double lambda, double nu, bool take_log, double hybrid_tol,
	double truncate_tol, double ymax);

//' Compute normalizing constant by truncating the infinite series.
//' 
//' @param lambda A vector of lambda parameters.
//' @param nu A vector of nu parameters.
//' @param tol Tolerance for truncation.
//' @param take_log Compute value on the log-scale if \code{TRUE}. Otherwise,
//' compute value on the original scale.
//' @param ymax Maximum value of y to consider.
//' 
//' @return Vector of normalizing constant values.
//' 
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector z_trunc(const Rcpp::NumericVector& lambda,
	const Rcpp::NumericVector& nu, double tol, bool take_log,
	double ymax);

//' An approximation method to compute normalizing constant (Shmueli et al,
//' JRSS-C, 2005).
//' 
//' @param lambda A vector of lambda parameters.
//' @param nu A vector of nu parameters.
//' @param take_log Compute value on the log-scale if \code{TRUE}. Otherwise,
//' compute value on the original scale.
//' 
//' @return Vector of normalizing constant values.
//' @references
//' Galit Shmueli, Thomas P. Minka, Joseph B. Kadane, Sharad Borle, and Peter
//' Boatwright (2005). A useful distribution for fitting discrete data: revival
//' of the Conway-Maxwell-Poisson distribution. Journal of Royal Statistical
//' Society C, 54, 127-142.
//' 
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector z_approx(const Rcpp::NumericVector& lambda,
	const Rcpp::NumericVector& nu, bool take_log);

//' Hybrid method to compute normalizing constant.
//' 
//' @param lambda A vector of lambda parameters.
//' @param nu A vector of nu parameters.
//' @param take_log Return constant on the log-scale if \code{TRUE}. Otherwise,
//' return constant on the original scale.
//' @param hybrid_tol Tolerance for when to use approximation vs. truncation.
//' @param truncate_tol Tolerance for truncation.
//' @param ymax Maximum value of y to consider.
//' 
//' @details
//' Conway-Maxwell Poisson normalizing constant is
//' \deqn{
//' Z(\lambda, \nu) = \sum_{j=0}^1 { \lambda^j / (j!)^\nu }.
//' }
//' Hybrid method uses approximation approach when
//' \deqn{
//' \lambda^{-1/\nu}
//' }
//' is smaller than \code{hybrid_tol}.
//' 
//' @return Vector of normalizing constant values.
//' 
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector z_hybrid(const Rcpp::NumericVector& lambda,
	const Rcpp::NumericVector& nu, bool take_log, double hybrid_tol,
	double truncate_tol, double ymax);

//' Truncation value computed in \code{z_trunc}. This can be used by other
//' functions to compute truncated expressions.
//' 
//' @param lambda A vector of lambda parameters.
//' @param nu A vector of nu parameters.
//' @param tol Tolerance for truncation.
//' @param ymax Maximum value of y to consider.
//' 
//' @return Vector of truncation values used in normalizing constant
//' computation.
//' 
//' @noRd
// [[Rcpp::export]]
Rcpp::IntegerVector y_trunc(const Rcpp::NumericVector& lambda,
	const Rcpp::NumericVector& nu, double tol, double ymax);

#endif
