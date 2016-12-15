glm.cmp <- function(formula, initial.est=NULL, nuinit=1, max=100, ...){
	initial_glm = glm(formula, family='poisson', ...)
	if (is.null(initial.est)) {
		initial.est = coef(initial_glm)
	}

	# Parse formula
	mf <- model.frame(formula, ...)
	y <- model.response(mf)
	X <- model.matrix(formula, mf)

	object_result = list()
	object_result$formula = formula
	object_result$response = y
	object_result$predictors = X
	object_result$max = max

	internal_result = ComputeBetasAndNuHat(object_result$predictors, object_result$response, betainit=initial.est, nuinit=nuinit, max=object_result$max)
	if(internal_result$convergence == 1) {stop(sprintf("Constant CMP estimates could not be determined.  Optimization scheme did not converge: ", internal_result))}

	num_pars = length(internal_result$par)
	object_result$glm_coefficients = coef(initial_glm)
	object_result$coefficients = internal_result$par[1:(num_pars -1)]
	object_result$nu = internal_result$par[num_pars]
	object_result$loglik <- -internal_result$objective
	object_result$convergence <- internal_result$convergence
	object_result$message <- internal_result$message

	attr(object_result, "class") <- c("cmp", attr(object_result, "class"))
	return(object_result)
}

coef.cmp <- function(object, ...) {
	return(object$coefficients)
}

nu.cmp <- function(object, ...) {
	return(object$nu)
}

sdev.cmp <- function(object, ...) {
	object$sdev <- CMPStdErrors(object$predictors, object$coef, object$nu, max=object$max)
	names(object$sdev) = c(colnames(object$predictors), "nu")
	return(object$sdev)
}

equitest.cmp <- function(object, ...) {
	res <- LRT.cmp(object$predictors, object$response, object$glm_coefficients, object$coef, object$nu, object$max)
	list(teststat = res$teststat[1,1], pvalue = res$pvalue)
}

leverage.cmp <- function(object, ...) {
	CMPLeverage(object$predictors, object$response, object$coef, object$nu, object$max)
}

deviance.cmp <- function(object, ...) {
	CMPDeviance(object$predictors, object$response, object$coef, object$nu, leverage.cmp(object), object$max)
}

residuals.cmp <- function(object, type = c("raw", "quantile"), ...) {
	X <- object$predictors
	lambda.hat <- exp(X %*% coef(object))
	nu.hat <- nu(object)
	y.hat <- predict(object, newdata = object$predictors)

	type <- match.arg(type)
	if (type == "raw") {
		res <- object$response - y.hat
	} else if (type == "quantile") {
		res <- rqres.cmp(object$response, lambda = lambda.hat, nu = nu.hat, max = object$max)
	} else {
		stop("Unsupported residual type")
	}

	return(as.numeric(res))
}

predict.cmp <- function(object, newdata = NULL, ...)
{
	if (!is.null(newdata)) {
		# If any of the original models had an intercept added via model.matrix, they
		# will have an "(Intercept)" column. Let's add an "(Intercept)" to newdata
		# in case the user didn't make one.
		newdata <- as.data.frame(newdata)
        newdata$'(Intercept)' <- 1
		X <- as.matrix(newdata[,colnames(object$predictors)])
	} else {
		X <- object$predictors
	}

	out <- constantCMPfitsandresids(object$coef, object$nu, X)
	y.hat <- out$fit
	return(y.hat)
}

parametric_bootstrap.cmp <- function(object, reps = 1000, report.period = reps+1, ...) {
	bootstrap_results = as.data.frame(CMPParamBoot(x=object$predictors, object$glm_coefficients, betahat=object$coef, nuhat=object$nu, n=reps, report.period=report.period, max=object$max)$CMPresult)
	names(bootstrap_results) = c(colnames(object$predictors), "nu", recursive=TRUE)
	return(bootstrap_results)
}
