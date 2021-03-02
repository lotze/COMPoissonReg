#include "util.h"
#include <vector>

// Solaris gives errors on CRAN if we do not define this.
double log(unsigned int x)
{
	return log(double(x));
}

// A version of the "which" function, since it does not appear to be
// provided in Rcpp Sugar yet.
Rcpp::IntegerVector which(const Rcpp::LogicalVector& x)
{
	std::vector<int> idx;
	for (int i = 0; i < x.size(); i++) {
		if (x(i)) {
			idx.push_back(i);
		}
	}

	return Rcpp::IntegerVector(idx.begin(), idx.end());
}

// This is intended to be like the linspace function in Armadillo
Rcpp::NumericVector linspace(double start, double end, unsigned int N)
{
	Rcpp::NumericVector out(N);
	double inc = (end - start) / (N-1);

	for (unsigned int i = 0; i < N; i++) {
		out(i) = start + i * inc;
	}

	return out;
}

// Compute quantiles of a discrete distribution with values 0, 1, ..., k-1 and
// associated *cumulative* probabilities cp(0), cp(1), ..., cp(k-1). Use a
// bisection search because cp can be rather large. q and cp can be given on
// the log-scale or probability scale, but they are expected to be compatible.
unsigned int qdiscrete(double q, const Rcpp::NumericVector& cp)
{
	unsigned int k = cp.size();
	if (q > cp(k-1)) {
		Rcpp::stop("q > max(cp)");
	}

	// Otherwise do a binary search
	unsigned int x_lo = 0;
	unsigned int x_hi = k-1;
	unsigned int x = int(floor((x_hi + x_lo) / 2.0));
	while (x_hi - x_lo > 1) {
		bool ind = (cp(x) >= q);
		x_lo = (1 - ind) * x + ind * x_lo;
		x_hi = ind * x + (1 - ind) * x_hi;
		x = int(floor((x_hi + x_lo) / 2.0));
	}

	if (cp(x_lo) >= q) {
		return x_lo;
	} else {
		return x_hi;
	}
}
