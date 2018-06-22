#include <Rcpp.h>
#include <list>

// A version of the "which" function, since it does not appear to be
// provided in Rcpp Sugar yet.
Rcpp::IntegerVector which(const Rcpp::LogicalVector& x)
{
	std::list<int> idx;
	for (int i = 0; i < x.size(); i++) {
		if (x(i)) {
			idx.push_back(i);
		}
	}

	return Rcpp::IntegerVector(idx.begin(), idx.end());
}

// Enumerate the terms lambda^y / (y!)^nu for y >= 0, until they become small,
// or until y = ymax is reached.
// [[Rcpp::export]]
Rcpp::NumericVector cmp_allprobs(double lambda, double nu, double tol = 1e-6,
	bool take_log = false, double ymax = 1e100, bool normalize = true)
{
	double z = 0;
	std::list<double> logp_unnorm;

	double delta;
	double psi = -0.5772157;	// Euler–Mascheroni constant
	double y = 0;
	double deriv = R_PosInf;
	while (deriv > 0 && y <= ymax && !isinf(z)) {
		delta = y*log(lambda) - nu*lgamma(y+1);
		logp_unnorm.push_back(delta);
		z += exp(delta);
		psi += 1 / (y+1);
		deriv = log(lambda) - nu*psi;
		y++;
	}

	double log_tol = log(tol);
	while (delta > log_tol && y <= ymax && !isinf(z)) {
		delta = y*log(lambda) - nu*lgamma(y+1);
		logp_unnorm.push_back(delta);
		z += exp(delta);
		y++;
	}

	Rcpp::NumericVector logp(logp_unnorm.begin(), logp_unnorm.end());
	if (normalize) {
		logp = logp - log(z);
	}

	if (take_log) {
		return logp;
	} else {
		return exp(logp);
	}
}

// Brute force computation of the z-function. This computation is done at the
// original scale (as opposed to the log-scale), so it becomes unstable when
// the magnitudes of the terms become very large.
// [[Rcpp::export]]
Rcpp::NumericVector z_exact(const Rcpp::NumericVector& lambda,
	const Rcpp::NumericVector& nu, double tol = 1e-6, bool take_log = false,
	double ymax = 1e100)
{
	unsigned int n = lambda.size();
	if (n != nu.size()) {
		Rcpp::stop("Length of lambda must be same as length of nu");
	}

	Rcpp::NumericVector z(n);
	for (unsigned int i = 0; i < n; i++) {
		const Rcpp::NumericVector& probs = cmp_allprobs(lambda(i), nu(i), tol, false, ymax, false);
		z(i) = Rcpp::sum(probs);
	}

	if (take_log) {
		return log(z);
	} else {
		return z;
	}
}

// Approximation to the z-function from Shmueli et al (JRSS-C, 2005)
// [[Rcpp::export]]
Rcpp::NumericVector z_approx(const Rcpp::NumericVector& lambda,
	const Rcpp::NumericVector& nu, bool take_log = false)
{
	unsigned int n = lambda.size();
	if (n != nu.size()) {
		Rcpp::stop("Length of lambda must be same as length of nu");
	}

	Rcpp::NumericVector logz = nu*exp(1/nu * log(lambda)) -
		(nu-1)/(2*nu) * log(lambda) -
		(nu-1)/2 * log(2*M_PI) - 0.5*log(nu);

	if (take_log) {
		return logz;
	} else {
		return exp(logz);
	}
}

// Sum (j>=0) { lambda^j / (j!)^nu }
// If lambda^{-1/nu} is small, it is more efficient (and often numerically more stable)
// to use an approximation at the log-scale.
// [[Rcpp::export]]
Rcpp::NumericVector z_hybrid(const Rcpp::NumericVector& lambda,
	const Rcpp::NumericVector& nu, bool take_log = false,
	double tol1 = 1e-2, double tol2 = 1e-6)
{
	unsigned int n = lambda.size();
	if (n != nu.size()) {
		Rcpp::stop("lambda and nu must be the same length");
	}

	Rcpp::LogicalVector cond = exp(-1/nu * log(lambda)) < tol1;
	Rcpp::IntegerVector idx_approx = which(cond);
	Rcpp::IntegerVector idx_exact = which(!cond);
	// Rcpp::IntegerVector idx_exact = Rcpp::setdiff(Rcpp::seq_len(n)-1, idx_approx);

	Rcpp::NumericVector logz(n);
	logz[idx_approx] = z_approx(lambda[idx_approx], nu[idx_approx], true);
	logz[idx_exact] = z_exact(lambda[idx_exact], nu[idx_exact], tol2, true);

	if (take_log) {
		return logz;
	} else {
		return exp(logz);
	}
}

// [[Rcpp::export]]
Rcpp::NumericVector cmp_expected_value(const Rcpp::NumericVector& lambda,
	const Rcpp::NumericVector& nu, double tol = 1e-6)
{
	unsigned int n = lambda.size();
	if (n != nu.size()) {
		Rcpp::stop("lambda and nu must be the same length");
	}
	
	Rcpp::NumericVector ex(n);
	ex.fill(0);

	for (unsigned int i = 0; i < n; i++) {
		Rcpp::NumericVector allprobs = cmp_allprobs(lambda(i), nu(i), tol);
		// Rcpp::IntegerVector values = Rcpp::seq_len(allprobs.size()) - 1;
		// ex(i) = sum(allprobs * values);

		for (unsigned int j = 0; j < allprobs.size(); j++) {
			ex(i) += j * allprobs(j);
		}
	}

	return ex;
}

// Sum (from j=0 to j=max) of j*lambda^j/((j!)^nu)
/*
// [[Rcpp::export]]
Rcpp::NumericVector computez_prodj_cpp(const Rcpp::NumericVector& lambda,
	const Rcpp::NumericVector& nu, unsigned int max = 100)
{
	unsigned int n = lambda.size();
	if (n != nu.size()) {
		Rcpp::stop("lambda and nu must be the same length");
	}

	Rcpp::NumericVector res(n);
	res.fill(0.0);

	for (unsigned int j = 0; j < max; j++) {
		res += exp(log(j) + j*log(lambda) - nu*lgamma(j+1));
	}

	return res;
}
*/

// [[Rcpp::export]]
Rcpp::NumericVector dcmp_cpp(const Rcpp::NumericVector& x,
	const Rcpp::NumericVector& lambda, const Rcpp::NumericVector& nu,
	double tol = 1e-6, bool take_log = false)
{
	unsigned int n = x.size();
	if (n != lambda.size() || n != nu.size()) {
		Rcpp::stop("lambda and nu must have length n");
	}

	Rcpp::NumericVector fx(n);
	fx.fill(0);

	for (unsigned int i = 0; i < n; i++) {
		Rcpp::NumericVector allprobs = cmp_allprobs(lambda(i), nu(i), tol);

		// Rcpp::print(allprobs);
		// Rprintf("Looking for element with index %d", x(i) + 1);

		if (x(i) <= allprobs.size() - 1) {
			fx(i) = allprobs(x(i));
		}
	}

	if (take_log) {
		return log(fx);
	} else {
		return fx;
	}
}

// [[Rcpp::export]]
Rcpp::NumericVector pcmp_cpp(const Rcpp::NumericVector& x,
	const Rcpp::NumericVector& lambda, const Rcpp::NumericVector& nu,
	double tol = 1e-6)
{
	unsigned int n = x.size();
	if (n != lambda.size() || n != nu.size()) {
		Rcpp::stop("lambda and nu must have length n");
	}

	Rcpp::NumericVector Fx(n);
	Fx.fill(1);

	for (unsigned int i = 0; i < n; i++) {
		Rcpp::NumericVector allprobs = cmp_allprobs(lambda(i), nu(i), tol);
		Rcpp::NumericVector cumprobs = Rcpp::cumsum(allprobs);
		if (x(i) <= allprobs.size() - 1) {
			Fx(i) = cumprobs(x(i));
		}
	}

	return Fx;
}

// Assume that quantiles q are not given on the log scale
// [[Rcpp::export]]
Rcpp::NumericVector qcmp_cpp(const Rcpp::NumericVector& q,
	const Rcpp::NumericVector& lambda, const Rcpp::NumericVector& nu,
	double tol = 1e-6)
{
	unsigned int n = q.size();
	if (n != lambda.size() || n != nu.size()) {
		Rcpp::stop("lambda and nu must have length n");
	}

	Rcpp::NumericVector x(n);

	for (unsigned int i = 0; i < n; i++) {
		Rcpp::NumericVector allprobs = cmp_allprobs(lambda(i), nu(i), tol);
		Rcpp::NumericVector cumprobs = Rcpp::cumsum(allprobs);
		Rcpp::IntegerVector idx = which(q(i) <= cumprobs);

		// Rcpp::print(cumprobs);
		// Rprintf("Looking for quantile %f", q(i));
		// Rcpp::print(idx);

		if (idx.size() > 0) {
			x(i) = idx(0);
		} else {
			x(i) = cumprobs.size();
		}
	}

	return x;
}

// Produce n(i) iid draws for each lambda(i) and nu(i).
// If n is a scalar, take n(i) to be n
// [[Rcpp::export]]
Rcpp::NumericVector rcmp_cpp(unsigned int n, const Rcpp::NumericVector& lambda,
	const Rcpp::NumericVector& nu, double tol = 1e-6)
{
	if (n != lambda.size() || n != nu.size()) {
		Rcpp::stop("lambda and nu must have length n");
	}

	Rcpp::NumericVector u = Rcpp::runif(n, 0.0, 1.0);
	return qcmp_cpp(u, lambda, nu, tol);
}
