glm.cmp <- function(formula.lambda, formula.nu = ~ 1, formula.p = NULL,
	beta.init = NULL, gamma.init = NULL, zeta.init = NULL, max = 100, ...)
{
	# Parse formula.lambda. This one should have the response.
	mf <- model.frame(formula.lambda, ...)
	y <- model.response(mf)
	X <- model.matrix(formula.lambda, mf)
	d1 <- ncol(X)

	# Parse formula.nu
	if (is.null(formula.nu)) {
		stop("formula.nu must be specified (can not be NULL)")
	}
	mf <- model.frame(formula.nu, ...)
	S <- model.matrix(formula.nu, mf)
	d2 <- ncol(S)

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
	res$max <- max
	res$beta.init <- beta.init
	res$gamma.init <- gamma.init

	# Handle ZI and non-ZI cases separately.
	if (!is.null(formula.p)) {
		mf <- model.frame(formula.p, ...)
		W <- model.matrix(formula.p, mf)
		d3 <- ncol(W)
		res$W <- W

		if (is.null(zeta.init)) { zeta.init <- rep(0, d3) }

		fit.out <- fit.zicmp.reg(res$y, res$X, res$S, res$W, beta.init = beta.init,
			gamma.init = gamma.init, zeta.init = zeta.init, max = res$max)

		res$zeta.init <- zeta.init
		res$beta.glm <- coef(initial.glm)
		res$beta <- fit.out$theta.hat$beta
		res$gamma <- fit.out$theta.hat$gamma
		res$zeta <- fit.out$theta.hat$zeta
		res$FIM <- fit.out$FIM
		res$V <- fit.out$V
		res$loglik <- fit.out$loglik
		res$opt.res <- fit.out$opt.res
		res$opt.method <- getOption("COMPoissonReg.optim.method")
		res$elapsed.sec <- fit.out$elapsed.sec

		attr(res, "class") <- c("zicmp", attr(res, "class"))
	} else {
		fit.out <- fit.cmp.reg(res$y, res$X, res$S, beta.init = beta.init,
			gamma.init = gamma.init, max = res$max)

		res$beta.glm <- coef(initial.glm)
		res$beta <- fit.out$theta.hat$beta
		res$gamma <- fit.out$theta.hat$gamma
		res$FIM <- fit.out$FIM
		res$V <- fit.out$V
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