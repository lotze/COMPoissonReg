fitted_zip_internal = function(X, W, beta, zeta, off_x, off_w)
{
	list(
		lambda = as.numeric(exp(X %*% beta + off_x)),
		p = as.numeric(plogis(W %*% zeta + off_w))
	)
}
