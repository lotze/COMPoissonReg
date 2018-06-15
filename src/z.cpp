#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::NumericVector computez_adapt(const Rcpp::NumericVector& lambda,
	const Rcpp::NumericVector& nu, double tol = 1e-6, bool take_log = false,
	double ymax = 1e100)
{
	unsigned int n = lambda.size();
	if (n != nu.size()) {
		Rcpp::stop("Length of lambda must be same as length of nu");
	}
	Rcpp::NumericVector z(n);
	z.fill(0);

	for (unsigned int i = 0; i < n; i++) {
		double delta = R_PosInf;
		double psi = -0.5772157;	// Euler–Mascheroni constant
		double y = 0;
		double deriv = R_PosInf;
		while (deriv > 0 && y <= ymax && !isinf(z(i))) {
			delta = exp( y*log(lambda(i)) - nu(i)*lgamma(y+1) );
			z(i) += delta;
			psi += 1 / (y+1);
			deriv = log(lambda(i)) - nu(i)*psi;
			y++;
		}

		while (delta > tol && y <= ymax && !isinf(z(i))) {
			delta = exp( y*log(lambda(i)) - nu(i)*lgamma(y+1) );
			z(i) += delta;
			y++;
		}
	}

	if (take_log) {
		return log(z);
	} else {
		return z;
	}
}
