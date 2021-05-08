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
* Compute quantiles of a discrete distribution with values 0, 1, ..., k-1 and
* associated *cumulative* probabilities cp(0), cp(1), ..., cp(k-1). Use a
* bisection search in case cp is a large vector . q and cp can be given on
* the log-scale or probability scale, but they are expected to be compatible.
*/
// [[Rcpp::export]]
unsigned int q_discrete(double q, const Rcpp::NumericVector& cp);

// Solaris gave errors on CRAN if we do not define this.
double log(unsigned int x);

// A version of R `which` function for Rcpp
Rcpp::IntegerVector which(const Rcpp::LogicalVector& x);

#endif
