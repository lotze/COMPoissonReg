glm.cmp <- function(formula.lambda, formula.nu = NULL, formula.p = NULL,
	beta.init = NULL, gamma.init = NULL, zeta.init = NULL, ...)
{
	# Parse formula.lambda. This one should have the response.
	mf <- model.frame(formula.lambda, ...)
	y <- model.response(mf)
	X <- model.matrix(formula.lambda, mf)
	off.X <- model.offset(mf)
	d1 <- ncol(X)

	# Parse formula.nu
	if (is.null(formula.nu)) {
		formula.nu = y ~ 1
	}
	mf <- model.frame(formula.nu, ...)
	S <- model.matrix(formula.nu, mf)
	off.S <- model.offset(mf)
	d2 <- ncol(S)

	if (is.null(off.X)) { off.X <- rep(0, n) }
	if (is.null(off.S)) { off.S <- rep(0, n) }

	initial.glm <- glm(formula.lambda, family='poisson', ...)
	if (is.null(beta.init)) { beta.init <- coef(initial.glm) }
	if (is.null(gamma.init)) { gamma.init <- rep(0, d2) }

	res <- list()
	res$formula.lambda <- formula.lambda
	res$formula.nu <- formula.nu
	res$formula.p <- formula.p
	res$y <- y
	res$X <- X
	res$S <- S
	res$beta.init <- beta.init
	res$gamma.init <- gamma.init
	res$off.X <- off.X
	res$off.S <- off.S

	# Handle ZI and non-ZI cases separately.
	if (!is.null(formula.p)) {
		mf <- model.frame(formula.p, ...)
		W <- model.matrix(formula.p, mf)
		off.W <- model.offset(mf)
		if (is.null(off.W)) { off.W <- rep(0, n) }
		d3 <- ncol(W)
		res$W <- W
		res$off.W <- off.W

		if (is.null(zeta.init)) { zeta.init <- rep(0, d3) }

		fit.out <- fit.zicmp.reg(res$y, res$X, res$S, res$W,
			beta.init = beta.init, gamma.init = gamma.init, zeta.init = zeta.init,
			off.X = off.X, off.S = off.S, off.W = off.W)

		res$zeta.init <- zeta.init
		res$beta.glm <- coef(initial.glm)
		res$beta <- fit.out$theta.hat$beta
		res$gamma <- fit.out$theta.hat$gamma
		res$zeta <- fit.out$theta.hat$zeta
		res$H <- fit.out$H
		res$loglik <- fit.out$loglik
		res$opt.res <- fit.out$opt.res
		res$opt.method <- getOption("COMPoissonReg.optim.method")
		res$elapsed.sec <- fit.out$elapsed.sec

		attr(res, "class") <- c("zicmp", attr(res, "class"))
	} else {
		fit.out <- fit.cmp.reg(res$y, res$X, res$S, beta.init = beta.init,
			gamma.init = gamma.init, off.X = off.X, off.S = off.S)

		res$beta.glm <- coef(initial.glm)
		res$beta <- fit.out$theta.hat$beta
		res$gamma <- fit.out$theta.hat$gamma
		res$H <- fit.out$H
		res$loglik <- fit.out$loglik
		res$opt.res <- fit.out$opt.res
		res$opt.method <- getOption("COMPoissonReg.optim.method")
		res$elapsed.sec <- fit.out$elapsed.sec

		attr(res, "class") <- c("cmp", attr(res, "class"))
	}

	# Add the test for equidispersion
	res$equitest <- equitest(res)

	return(res)
}