fitted_zip_internal = function(X, W, beta, zeta, off.x, off.w)
{
	list(
		lambda = as.numeric(exp(X %*% beta + off.x)),
		p = as.numeric(plogis(W %*% zeta + off.w))
	)
}
