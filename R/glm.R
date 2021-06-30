#' COM-Poisson and Zero-Inflated COM-Poisson Regression
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
#' @param off.x Offset term for the `lambda` regression. If \code{NULL}, offsets
#' are taken to be zeros.
#' @param off.s Offset term for the `nu` regression. If \code{NULL}, offsets
#' are taken to be zeros.
#' @param off.w Offset term for the `p` regression. If \code{NULL}, offsets
#' are taken to be zeros.
#' @param beta.init Initial value for `beta`. If \code{NULL}, initial value
#' is taken to be a vector of zeros.
#' @param gamma.init Initial value for `gamma`. If \code{NULL}, initial value
#' is taken to be a vector of zeros.
#' @param zeta.init Initial value for `zeta`. If \code{NULL}, initial value
#' is taken to be a vector of zeros.
#' 
#' @return
#' See the \link{glm.cmp}.
#' 
#' @name glm.cmp-raw
NULL

#' @name glm.cmp-raw
#' @export
get.fixed = function(beta = integer(0), gamma = integer(0), zeta = integer(0))
{
	stopifnot(is.numeric(beta))
	stopifnot(is.numeric(gamma))
	stopifnot(is.numeric(zeta))

	out = list(beta = beta, gamma = gamma, zeta = zeta)
	class(out) = "glm.cmp.fixed"
	return(out)
}

#' @name glm.cmp-raw
#' @export
get.init = function(beta = NULL, gamma = NULL, zeta = NULL)
{
	out = list(beta = beta, gamma = gamma, zeta = zeta)
	class(out) = "glm.cmp.init"
	return(out)
}

#' @name glm.cmp-raw
#' @export
get.offset = function(x, s, w)
{
	out = list(x = x, s = s, w = s)
	class(out) = "glm.cmp.offset"
	return(out)
}

#' @name glm.cmp
#' @export
glm.cmp = function(formula.lambda, formula.nu = ~ 1, formula.p = NULL,
	data = NULL, init = NULL, fixed = NULL, ...)
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
		off = get.offset(x = off.x, s = off.s)
		res = glm.cmp.raw(y, X, S, off = off, init = init, fixed = fixed)
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

		# Run the regression using the raw interface function
		off = get.offset(x = off.x, s = off.s, w = off.w)
		res = glm.zicmp.raw(y, X, S, W, off = off, init = init, fixed = fixed)
	}
	
	# Add a few things to return value
	res$interface = "formula"
	res$formula.lambda = formula.lambda
	res$formula.nu = formula.nu
	res$formula.p = formula.p

	return(res)
}

#' @name glm.cmp-raw
#' @export
glm.cmp.raw = function(y, X, S, off = NULL, init = NULL, fixed = NULL)
{
	# Get dimensions
	n = length(y)
	d1 = ncol(X)
	d2 = ncol(S)

	# Offsets
	if (is.null(off)) {
		off = get.offset(x = numeric(n), s = numeric(n), w = NULL)
	}
	stopifnot(class(off) == "glm.cmp.offset")

	# Fixed values
	if (is.null(fixed)) {
		fixed = get.fixed()
	}
	stopifnot(class(fixed) == "glm.cmp.fixed")

	# Make sure fixed indices are between 1 and the corresponding dimension
	stopifnot(all(fixed$beta %in% seq_len(d1)))
	stopifnot(all(fixed$gamma %in% seq_len(d2)))

	# Make sure dimensions match up
	stopifnot(n == nrow(X))
	stopifnot(n == nrow(S))
	stopifnot(n == length(off$x))
	stopifnot(n == length(off$s))

	# Initial parameter values
	if (is.null(init)) { init = get.init() }
	if (is.null(init$beta)) { init$beta = numeric(d1) }
	if (is.null(init$gamma)) { init$gamma = numeric(d2) }
	stopifnot(class(init) == "glm.cmp.init")
	stopifnot(d1 == length(init$beta))
	stopifnot(d2 == length(init$gamma))

	# Fit the CMP regression model
	fit.out = fit.cmp.reg(y, X, S, beta.init = init$beta,
		gamma.init = init$gamma, off.x = off$x, off.s = off$s,
		fixed.beta = fixed$beta, fixed.gamma = fixed$gamma)

	# Construct return value
	res = list(
		y = y,
		X = X,
		S = S,
		beta.init = init$beta,
		gamma.init = init$gamma,
		off.x = off$x,
		off.s = off$s,
		beta = fit.out$theta.hat$beta,
		gamma = fit.out$theta.hat$gamma,
		H = fit.out$H,
		loglik = fit.out$loglik,
		opt.res = fit.out$opt.res,
		optim.method = fit.out$optim.method,
		optim.control = fit.out$optim.control,
		elapsed.sec = fit.out$elapsed.sec,
		fixed.beta = fit.out$fixed.beta,
		fixed.gamma = fit.out$fixed.gamma,
		unfixed.beta = fit.out$unfixed.beta,
		unfixed.gamma = fit.out$unfixed.gamma,
		interface = "raw"
	)
	attr(res, "class") = c("cmp", attr(res, "class"))

	# Add the equidispersion test
	res$equitest = equitest(res)

	return(res)
}

#' @name glm.cmp-raw
#' @export
glm.zicmp.raw = function(y, X, S, W, off = NULL, init = NULL, fixed = NULL)
{
	# Get dimensions
	n = length(y)
	d1 = ncol(X)
	d2 = ncol(S)
	d3 = ncol(W)

	# Offsets
	if (is.null(off)) {
		off = get.offset(x = numeric(n), s = numeric(n), w = numeric(n))
	}
	stopifnot(class(off) == "glm.cmp.offset")

	# Fixed values
	if (is.null(fixed)) {
		fixed = get.fixed()
	}
	stopifnot(class(fixed) == "glm.cmp.fixed")

	# Make sure fixed indices are between 1 and the corresponding dimension
	stopifnot(all(fixed$beta %in% seq_len(d1)))
	stopifnot(all(fixed$gamma %in% seq_len(d2)))
	stopifnot(all(fixed$zeta %in% seq_len(d3)))

	# Make sure dimensions match up
	stopifnot(n == nrow(X))
	stopifnot(n == nrow(S))
	stopifnot(n == nrow(W))
	stopifnot(n == length(off$x))
	stopifnot(n == length(off$s))
	stopifnot(n == length(off$w))

	# Initial parameter values
	if (is.null(init)) { init = get.init() }
	if (is.null(init$beta)) { init$beta = numeric(d1) }
	if (is.null(init$gamma)) { init$gamma = numeric(d2) }
	if (is.null(init$zeta)) { init$zeta = numeric(d3) }

	stopifnot(class(init) == "glm.cmp.init")
	stopifnot(d1 == length(init$beta))
	stopifnot(d2 == length(init$gamma))
	stopifnot(d3 == length(init$zeta))

	# Fit the ZICMP regression model
	fit.out = fit.zicmp.reg(y, X, S, W, beta.init = init$beta,
		gamma.init = init$gamma, zeta.init = init$zeta, off.x = off$x,
		off.s = off$s, off.w = off$w, fixed.beta = fixed$beta,
		fixed.gamma = fixed$gamma, fixed.zeta = fixed$zeta)

	# Construct return value
	res = list(
		y = y,
		X = X,
		S = S,
		W = W,
		beta.init = init$beta,
		gamma.init = init$gamma,
		zeta.init = init$zeta,
		off.x = off$x,
		off.s = off$s,
		off.w = off$w,
		beta = fit.out$theta.hat$beta,
		gamma = fit.out$theta.hat$gamma,
		zeta = fit.out$theta.hat$zeta,
		H = fit.out$H,
		loglik = fit.out$loglik,
		opt.res = fit.out$opt.res,
		optim.method = fit.out$optim.method,
		optim.control = fit.out$optim.control,
		elapsed.sec = fit.out$elapsed.sec,
		fixed.beta = fit.out$fixed.beta,
		fixed.gamma = fit.out$fixed.gamma,
		fixed.zeta = fit.out$fixed.zeta,
		unfixed.beta = fit.out$unfixed.beta,
		unfixed.gamma = fit.out$unfixed.gamma,
		unfixed.zeta = fit.out$unfixed.zeta,
		interface = "raw"
	)
	attr(res, "class") = c("zicmp", attr(res, "class"))

	# Add the equidispersion test
	res$equitest = equitest(res)

	return(res)
}
