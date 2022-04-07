#ifndef Z_DERIVS_H
#define Z_DERIVS_H

#include <Rcpp.h>

/*
* Sum (from j=0 to j=max) of j*lambda^j/((j!)^nu)
*/
Rcpp::NumericVector z_prodj(const Rcpp::NumericVector& lambda,
	const Rcpp::NumericVector& nu, unsigned int max);

/*
* Sum (from j=0 to j=max) of j^2*lambda^j/((j!)^nu)
*/
Rcpp::NumericVector z_prodj2(const Rcpp::NumericVector& lambda,
	const Rcpp::NumericVector& nu, unsigned int max);

/*
* Sum (from j=0 to j=max) of jlog(j!)*lambda^j/((j!)^nu)
*/
Rcpp::NumericVector z_prodjlogj(const Rcpp::NumericVector& lambda,
	const Rcpp::NumericVector& nu, unsigned int max);

//' Derivatives of the normalizing constant.
//' 
//' @param lambda A vector of lambda parameters.
//' @param nu A vector of nu parameters.
//' @param max Maximum value of y to consider.
//' 
//' @details
//' Derivatives of the normalizing constant computed by truncation.
//' \itemize{
//' \item \code{z_prodlogj} \eqn{\sum_{j=0}^{\text{max}} \log(j!) \lambda^j/((j!)^\nu)}
//' }
//' Several other functions are currently in the code but not used.
//' 
//' @return Vector of results.
//' 
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector z_prodlogj(const Rcpp::NumericVector& lambda,
	const Rcpp::NumericVector& nu, unsigned int max);

/*
* Sum (from j=0 to j=max) of (log(j!))^2 * lambda^j/((j!)^nu)
*/
Rcpp::NumericVector z_prodlogj2(const Rcpp::NumericVector& lambda,
	const Rcpp::NumericVector& nu, unsigned int max);

#endif
