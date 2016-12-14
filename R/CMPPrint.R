summary.cmp <- function(object, ...)
{
	est <- c(coef(object), nu(object))
	se <- sdev(object)
	z.val <- est / se
	p.val <- 2*(1 - pnorm(abs(est / se)))

	DF <- data.frame(
		Estimate = round(est, 4),
		SE = round(se, 4),
		z.value = round(z.val, 6),
		p.value = sprintf("%0.4g", p.val)
	)
	rownames(DF) <- c(sprintf("X:%s", colnames(object$predictors)), "nu")

	list(DF = DF,
		n = nrow(object$predictors),
		loglik = logLik(object),
		aic = AIC(object),
		bic = BIC(object),
		opt.message = object$message,
		opt.convergence = object$convergence,
		elapsed.sec = object$elapsed.sec
	)
}

print.cmp <- function(x, ...)
{
	cat("Fit for CMP model\n")
	s <- summary(x)
	print(s$DF)

	cat("--\n")
	cat(sprintf("Elapsed Sec: %0.2f   ", s$elapsed.sec))
	cat(sprintf("Sample size: %d\n", s$n))
	cat(sprintf("LogLik: %0.4f   ", s$loglik))
	cat(sprintf("AIC: %0.4f   ", s$aic))
	cat(sprintf("BIC: %0.4f   ", s$bic))
	cat("\n")
	cat(sprintf("Converged status: %d   ", s$opt.convergence))
	cat(sprintf("Message: %s\n", s$opt.message))
}

logLik.cmp <- function(object, ...)
{
	object$loglik
}

AIC.cmp <- function(object, k, ...)
{
	-2*object$loglik + 2*length(coef(object))
}

BIC.cmp <- function(object, ...)
{
	n <- length(object$response)
	-2*object$loglik + log(n)*length(coef(object))
}

