#ifndef UTIL_H
#define UTIL_H

#include <Rcpp.h>

/*
* The following macro helps us to use the property:
*  log(x + y) = log(x) + log(1 + y/x)
*             = log(x) + log1p(exp(log(y) - log(x))),
*  with log(x) and log(y) given as inputs, i.e. on the log-scale. When x and y
*  are of very different magnitudes, this is more stable when x is taken to be
*  the larger of the inputs. The most extreme case is when one of the inputs
*  might be -inf; in this case that input should be the second one.
*  
*  https://en.wikipedia.org/wiki/List_of_logarithmic_identities#Summation
*/
#define logadd(logx, logy) (logx + log1p(exp(logy - logx)))

/*
* Solaris gave errors on CRAN if we did not define this.
*/
double log(unsigned int x);

/*
* A version of R `which` function for Rcpp.
*/
Rcpp::IntegerVector which(const Rcpp::LogicalVector& x);

//' Compute quantiles of a discrete / finite distribution.
//' 
//' @param q Probability of the desired quantile.
//' @param cp vector of cumulative probabilities.
//' 
//' @details
//' Compute a quantile for the discrete distribution with values
//' \code{0, 1, ..., k-1} and associated *cumulative* probabilities
//' \code{cp(0), cp(1), ..., cp(k-1)}. Use a bisection search in case \code{cp}
//' is a large vector. \code{q} and \code{cp} can be given on the log-scale or
//' probability scale, but they are expected to be compatible.
//' 
//' @return The desired quantile.
//' 
//' @noRd
// [[Rcpp::export]]
unsigned int q_discrete(double q, const Rcpp::NumericVector& cp);

#endif
