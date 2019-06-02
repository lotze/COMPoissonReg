fitted_zip_internal = function(X, W, beta, zeta, off.X, off.W)
{
	list(
		lambda = as.numeric(exp(X %*% beta + off.X)),
		p = as.numeric(plogis(W %*% zeta + off.W))
	)
}
