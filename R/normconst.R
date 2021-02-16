#' TBD: Update this description
#' TBD: Pass around options instead of relying on defaults.
#' Brute force computation of the z-function. This computation is done at the
#' original scale (as opposed to the log-scale), so it becomes unstable when
#' the magnitudes of the terms become very large.
#' Sum (j>=0) { lambda^j / (j!)^nu }
#' @name CMP Normalizing Constant
#' @export
ncmp = function(lambda, nu, log = FALSE, method = "hybrid")
{
	ymax = getOption("COMPoissonReg.ymax")
	hybrid_tol = getOption("COMPoissonReg.hybrid_tol")
	truncate_tol = getOption("COMPoissonReg.truncate_tol")

	n = max(length(lambda), length(nu))
	prep = prep.zicmp(n, lambda, nu)

	if (method == "trunc") {
		out = z_trunc(prep$lambda, prep$nu, tol = truncate_tol, take_log = log,
			ymax = ymax)
	} else if (method == "approx") {
		out = z_approx(prep$lambda, prep$nu, take_log = log)
	} else if (method == "hybrid") {
		out = z_hybrid(prep$lambda, prep$nu, take_log = log,
			hybrid_tol = hybrid_tol, truncate_tol = truncate_tol, ymax = ymax)
	} else {
		msg = sprintf("Method not recognized: %s", method)
		stop(msg)
	}

	return(out)
}
	
