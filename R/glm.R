#' COM-Poisson and Zero-Inflated COM-Poisson regression
#' 
#' Fit COM-Poisson regression using maximum likelihood estimation.
#' Zero-Inflated COM-Poisson can be fit by specifying a regression for the
#' overdispersion parameter.
#' 
#' @param formula.lambda regression formula linked to \code{log(lambda)}.
#'   The response should be specified here.
#' @param data An optional data.frame with variables to be used with regression
#'   formulas. Variables not found here are read from the envionment.
#' @param formula.nu regression formula linked to \code{log(nu)}. The
#'   default, is taken to be only an intercept.
#' @param formula.p regression formula linked to \code{logit(p)}. If NULL
#'   (the default), zero-inflation term is excluded from the model.
#' @param beta.init initial values for regression coefficients of \code{lambda}.
#' @param gamma.init initial values for regression coefficients of \code{nu}.
#' @param zeta.init initial values for regression coefficients of \code{p}.
#' @param ... other arguments, such as \code{subset} and \code{na.action}.
#' 
#' @return
#' \code{glm.cmp} produces an object of either class 'cmp' or 'zicmp', depending
#' on whether zero-inflation is used in the model. From this object, coefficients
#' and other information can be extracted.
#' 
#' @details 
#' The COM-Poisson regression model is
#' \deqn{
#' y_i \sim \rm{CMP}(\lambda_i, \nu_i), \;\;\;
#' \log \lambda_i = \bm{x}_i^\top \beta, \;\;\;
#' \log \nu_i = \bm{s}_i^\top \gamma.
#' }{
#' y_i ~ CMP(lambda_i, nu_i),
#' log lambda_i = x_i^T beta,
#' log nu_i = s_i^T gamma.
#' }
#' 
#' The Zero-Inflated COM-Poisson regression model assumes that \eqn{y_i} is 0
#' with probability \eqn{p_i} or \eqn{y_i^*} with probability \eqn{1 - p_i},
#' where
#' \deqn{
#' y_i^* \sim \rm{CMP}(\lambda_i, \nu_i), \;\;\;
#' \log \lambda_i = \bm{x}_i^\top \beta, \;\;\;
#' \log \nu_i = \bm{s}_i^\top \gamma, \;\;\;
#' \log p_i = \bm{w}_i^\top \zeta.
#' }{
#' y_i^* ~ CMP(lambda_i, nu_i),
#' log lambda_i = x_i^T beta,
#' log nu_i = s_i^T gamma,
#' log p_i = w_i^T zeta.
#' }
#'
#' @references
#' Kimberly F. Sellers & Galit Shmueli (2010). A Flexible Regression Model for
#' Count Data. Annals of Applied Statistics, 4(2), 943-961.
#' 
#' Kimberly F. Sellers and Andrew M. Raim (2016). A Flexible Zero-Inflated Model
#' to Address Data Dispersion. Computational Statistics and Data Analysis, 99,
#' 68-80.
#'
#' @author Kimberly Sellers, Thomas Lotze, Andrew Raim
#' @name glm.cmp
#' @export
glm.cmp = function(formula.lambda, formula.nu = ~ 1, formula.p = NULL,
	data = NULL, beta.init = NULL, gamma.init = NULL, zeta.init = NULL, ...)
{
	# Parse formula.lambda. This one should have the response.
	mf = model.frame(formula.lambda, data, ...)
	y = model.response(mf)
	X = model.matrix(formula.lambda, mf)
	off.x = model.offset(mf)
	d1 = ncol(X)

	weights = model.weights(mf)
	if(!is.null(weights)) {
		stop("weights argument is currently not supported")
	}

	# Parse formula.nu
	mf = model.frame(formula.nu, data, ...)
	S = model.matrix(formula.nu, mf)
	off.s = model.offset(mf)
	d2 = ncol(S)

	n = length(y)
	if (is.null(off.x)) { off.x = rep(0, n) }
	if (is.null(off.s)) { off.s = rep(0, n) }

	# Build up a call to glm function
	# Try to do it in a general way that mimics the call to glm.cmp;
	# e.g. data should come from the same environment, and ... arguments
	# like "subset" should be passed along.
	glm.args = list(
		formula = formula.lambda,
		family = poisson
	)
	if (!is.null(data)) { glm.args$data = data }
	glm.args = c(glm.args, list(...))
	initial.glm = do.call(stats::glm, glm.args)

	if (is.null(beta.init)) { beta.init = coef(initial.glm) }
	if (is.null(gamma.init)) { gamma.init = rep(0, d2) }

	res = list(
		formula.lambda = formula.lambda,
		formula.nu = formula.nu,
		formula.p = formula.p,
		y = y,
		X = X,
		S = S,
		beta.init = beta.init,
		gamma.init = gamma.init,
		off.x = off.x,
		off.s = off.s
	)

	# Handle ZI and non-ZI cases separately.
	if (!is.null(formula.p)) {
		mf = model.frame(formula.p, data, ...)
		W = model.matrix(formula.p, mf)
		off.w = model.offset(mf)
		if (is.null(off.w)) { off.w = rep(0, n) }
		d3 = ncol(W)
		res$W = W
		res$off.w = off.w

		if (is.null(zeta.init)) { zeta.init = rep(0, d3) }

		fit.out = fit.zicmp.reg(res$y, res$X, res$S, res$W,
			beta.init = beta.init, gamma.init = gamma.init, zeta.init = zeta.init,
			off.x = off.x, off.s = off.s, off.w = off.w)

		res$zeta.init = zeta.init
		res$beta.glm = coef(initial.glm)
		res$beta = fit.out$theta.hat$beta
		res$gamma = fit.out$theta.hat$gamma
		res$zeta = fit.out$theta.hat$zeta
		res$H = fit.out$H
		res$loglik = fit.out$loglik
		res$opt.res = fit.out$opt.res
		res$opt.method = getOption("COMPoissonReg.optim.method")
		res$elapsed.sec = fit.out$elapsed.sec

		attr(res, "class") = c("zicmp", attr(res, "class"))
	} else {
		fit.out = fit.cmp.reg(res$y, res$X, res$S, beta.init = beta.init,
			gamma.init = gamma.init, off.x = off.x, off.s = off.s)

		res$beta.glm = coef(initial.glm)
		res$beta = fit.out$theta.hat$beta
		res$gamma = fit.out$theta.hat$gamma
		res$H = fit.out$H
		res$loglik = fit.out$loglik
		res$opt.res = fit.out$opt.res
		res$opt.method = getOption("COMPoissonReg.optim.method")
		res$elapsed.sec = fit.out$elapsed.sec

		attr(res, "class") = c("cmp", attr(res, "class"))
	}

	# Add the test for equidispersion
	res$equitest = equitest(res)

	return(res)
}
