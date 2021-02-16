#ifndef Z_DERIVS_H
#define Z_DERIVS_H

#include <Rcpp.h>

// Sum (from j=0 to j=max) of j*lambda^j/((j!)^nu)
// [[Rcpp::export]]
Rcpp::NumericVector z_prodj(const Rcpp::NumericVector& lambda,
	const Rcpp::NumericVector& nu, unsigned int max = 100);

// Sum (from j=0 to j=max) of j^2*lambda^j/((j!)^nu)
// [[Rcpp::export]]
Rcpp::NumericVector z_prodj2(const Rcpp::NumericVector& lambda,
	const Rcpp::NumericVector& nu, unsigned int max = 100);

// Sum (from j=0 to j=max) of jlog(j!)*lambda^j/((j!)^nu)
// [[Rcpp::export]]
Rcpp::NumericVector z_prodjlogj(const Rcpp::NumericVector& lambda,
	const Rcpp::NumericVector& nu, unsigned int max = 100);

// Sum (from j=0 to j=max) of log(j!)*lambda^j/((j!)^nu)
// [[Rcpp::export]]
Rcpp::NumericVector z_prodlogj(const Rcpp::NumericVector& lambda,
	const Rcpp::NumericVector& nu, unsigned int max = 100);

// Sum (from j=0 to j=max) of (log(j!))^2 * lambda^j/((j!)^nu)
// [[Rcpp::export]]
Rcpp::NumericVector z_prodlogj2(const Rcpp::NumericVector& lambda,
	const Rcpp::NumericVector& nu, unsigned int max = 100);

#endif
