#' COM-Poisson and Zero-Inflated COM-Poisson Regression
#' 
#' Fit COM-Poisson regression using maximum likelihood estimation.
#' Zero-Inflated COM-Poisson can be fit by specifying a regression for the
#' overdispersion parameter.
#' 
#' @param formula.lambda regression formula linked to \code{log(lambda)}.
#'   The response should be specified here.
#' @param formula.nu regression formula linked to \code{log(nu)}. The
#'   default, is taken to be only an intercept.
#' @param formula.p regression formula linked to \code{logit(p)}. If NULL
#'   (the default), zero-inflation term is excluded from the model.
#' @param data An optional data.frame with variables to be used with regression
#'   formulas. Variables not found here are read from the envionment.
#' @param offset A data structure that specifies offsets. See the helper
#' function \link{get.offset}.
#' @param init A data structure that specifies initial values. See the helper
#' function \link{get.init}.
#' @param fixed A data structure that specifies which coefficients should
#' remain fixed in the maximum likelihood procedure. See the helper function
#' \link{get.fixed}.
#' @param control A control data structure. See the helper function
#' \link{get.control}. If \code{NULL}, a global default will be used.
#' @param ... other arguments, such as \code{subset} and \code{na.action}.
#' 
#' @return
#' \code{glm.cmp} produces an object of either class \code{cmpfit} or
#' \code{zicmpfit}, depending on whether zero-inflation is used in the model.
#' From this object, coefficients and other information can be extracted.
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
#' \rm{logit} \, p_i = \bm{w}_i^\top \zeta.
#' }{
#' y_i^* ~ CMP(lambda_i, nu_i),
#' log lambda_i = x_i^T beta,
#' log nu_i = s_i^T gamma,
#' logit p_i = w_i^T zeta.
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
NULL

#' Raw Interface to COM-Poisson and Zero-Inflated COM-Poisson Regression
#' 
#' Fit COM-Poisson and Zero-Inflated COM-Poisson regression using a "raw"
#' interface which bypasses the formula-driven interface of \code{glm.cmp}.
#' 
#' @param y A vector of counts which represent the response .
#' @param X Design matrix for the `lambda` regression.
#' @param S Design matrix for the `nu` regression.
#' @param W Design matrix for the `p` regression.
#' @param offset A data structure that specifies offsets. See the helper
#' function \link{get.offset}.
#' @param init A data structure that specifies initial values. See the helper
#' function \link{get.init}.
#' @param fixed A data structure that specifies which coefficients should
#' remain fixed in the maximum likelihood procedure. See the helper function
#' \link{get.fixed}.
#' @param control A control data structure. See the helper function
#' \link{get.control}.
#' 
#' @return
#' See the \link{glm.cmp}.
#' 
#' @name glm.cmp-raw
NULL

#' Extract model elements from a formula to use with the raw interface
#' 
#' @param formula.lambda regression formula linked to \code{log(lambda)}.
#' The response should be specified here.
#' @param formula.nu regression formula linked to \code{log(nu)}. The
#' default, is taken to be only an intercept.
#' @param formula.p regression formula linked to \code{logit(p)}. If NULL
#' (the default), zero-inflation term is excluded from the model.#'
#' @param data An optional data.frame with variables to be used with regression
#'   formulas. Variables not found here are read from the envionment.
#' @param ... other arguments, such as \code{subset} and \code{na.action}.
#'
#' @noRd
formula2raw = function(formula.lambda, formula.nu, formula.p, data = NULL, ...)
{
	# Parse formula.lambda. This one should have the response.
	mf = model.frame(formula.lambda, data, ...)
	y = model.response(mf)
	X = model.matrix(formula.lambda, mf)
	off.x = model.offset(mf)
	d1 = ncol(X)
	n = length(y)

	weights = model.weights(mf)
	if(!is.null(weights)) {
		stop("weights argument is currently not supported")
	}

	# Parse formula.nu
	mf = model.frame(formula.nu, data, ...)
	S = model.matrix(formula.nu, mf)
	if (nrow(S) == 0) {
		# A workaround for the case where there is no context in formula.nu
		# about how many observations there should be. The only way this
		# seems possible is when formula.nu = ~1
		S = model.matrix(y ~ 1)
	}
	off.s = model.offset(mf)
	d2 = ncol(S)

	if (is.null(off.x)) { off.x = rep(0, n) }
	if (is.null(off.s)) { off.s = rep(0, n) }

	# If formula.p is NULL, do CMP regression. Else do ZICMP regression
	if (is.null(formula.p)) {
		# Run the regression using the raw interface function
		W = matrix(0, n, 0)
		off.w = numeric(n)
	} else {
		mf = model.frame(formula.p, data, ...)
		W = model.matrix(formula.p, mf)
		if (nrow(W) == 0) {
			# A workaround for the case where there is no context in formula.nu
			# about how many observations there should be. The only way this
			# seems possible is when formula.p = ~1
			W = model.matrix(y ~ 1)
		}
		d3 = ncol(W)
		off.w = model.offset(mf)
		if (is.null(off.w)) { off.w = rep(0, n) }
	}

	offset = get.offset(x = off.x, s = off.s, w = off.w)
	res = list(y = y, X = X, S = S, W = W, offset = offset)
	return(res)
}

#' @name glm.cmp
#' @export
glm.cmp = function(formula.lambda, formula.nu = ~ 1, formula.p = NULL,
	data = NULL, init = NULL, fixed = NULL, control = NULL, ...)
{

	raw = formula2raw(formula.lambda, formula.nu, formula.p, data, ...)
	d3 = ncol(raw$W)

	if (d3 > 0) {
		res = glm.zicmp.raw(y = raw$y, X = raw$X, S = raw$S, W = raw$W,
			offset = raw$offset, init = init, fixed = fixed, control = control)
	} else {
		res = glm.cmp.raw(y = raw$y, X = raw$X, S = raw$S, offset = raw$offset,
			init = init, fixed = fixed, control = control)
	}

	# Add formula-specific things to return object
	res$interface = "formula"
	res$formula.lambda = formula.lambda
	res$formula.nu = formula.nu
	res$formula.p = formula.p

	return(res)
}

#' @name glm.cmp-raw
#' @export
glm.cmp.raw = function(y, X, S, offset = NULL, init = NULL, fixed = NULL, control = NULL)
{
	# Get dimensions
	n = length(y)
	d1 = ncol(X)
	d2 = ncol(S)

	# Initialize NULL arguments
	if (is.null(offset)) { offset = get.offset.zero(n) }
	if (is.null(init)) { init = get.init.zero(d1, d2) }
	if (is.null(fixed)) { fixed = get.fixed() }
	if (is.null(control)) { control = getOption("COMPoissonReg.control") }

	# Fit the CMP regression model
	fit.out = fit.cmp.reg(y, X, S, init = init, offset = offset, fixed = fixed,
		control = control)

	# Construct return value
	res = list(
		y = y,
		X = X,
		S = S,
		init = init,
		offset = offset,
		beta = fit.out$theta.hat$beta,
		gamma = fit.out$theta.hat$gamma,
		H = fit.out$H,
		loglik = fit.out$loglik,
		opt.res = fit.out$opt.res,
		control = fit.out$control,
		elapsed.sec = fit.out$elapsed.sec,
		fixed = fit.out$fixed,
		unfixed = fit.out$unfixed,
		interface = "raw"
	)
	attr(res, "class") = c("cmpfit", attr(res, "class"))

	# Add the equidispersion test
	res$equitest = equitest(res)

	return(res)
}

#' @name glm.cmp-raw
#' @export
glm.zicmp.raw = function(y, X, S, W, offset = NULL, init = NULL, fixed = NULL, control = NULL)
{
	# Get dimensions
	n = length(y)
	d1 = ncol(X)
	d2 = ncol(S)
	d3 = ncol(W)

	# Initialize NULL arguments
	if (is.null(offset)) { offset = get.offset.zero(n) }
	if (is.null(init)) { init = get.init.zero(d1, d2, d3) }
	if (is.null(fixed)) { fixed = get.fixed() }
	if (is.null(control)) { control = getOption("COMPoissonReg.control") }

	# Fit the ZICMP regression model
	fit.out = fit.zicmp.reg(y, X, S, W, init = init, offset = offset,
		fixed = fixed, control = control)

	# Construct return value
	res = list(
		y = y,
		X = X,
		S = S,
		W = W,
		init = init,
		offset = offset,
		beta = fit.out$theta.hat$beta,
		gamma = fit.out$theta.hat$gamma,
		zeta = fit.out$theta.hat$zeta,
		H = fit.out$H,
		loglik = fit.out$loglik,
		opt.res = fit.out$opt.res,
		control = fit.out$control,
		elapsed.sec = fit.out$elapsed.sec,
		fixed = fit.out$fixed,
		unfixed = fit.out$unfixed,
		interface = "raw"
	)
	attr(res, "class") = c("zicmpfit", attr(res, "class"))

	# Add the equidispersion test
	res$equitest = equitest(res)

	return(res)
}
