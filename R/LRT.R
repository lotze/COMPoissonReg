LRT.cmp <- function(x,y,betahat0,betahat,nuhat,max){
# This function computes the -2logLRT value and associated p-value for significance

# Compute the test statistic
  teststat <- -2*((t(y)%*%(x %*% betahat0)) - sum(lgamma(y+1)) - sum(exp(x %*% betahat0)) - 
              ((t(y)%*%(x %*% betahat)) - (nuhat*sum(lgamma(y+1))) - 
              sum(log(computez(exp(x %*% betahat),nuhat,max)))))

# Determine the associated p-value
  pvalue <- pchisq(teststat,df=1,lower.tail=FALSE)

return(list(teststat=teststat,pvalue=pvalue))
}

LRT.zicmp <- function(y, X, S, W, beta.hat, gamma.hat, zeta.hat, beta0.hat, zeta0.hat, max)
{
	n <- length(y)
	ff <- numeric(n)
	ff0 <- numeric(n)

	lambda.hat <- exp(X %*% beta.hat)
	nu.hat <- exp(S %*% gamma.hat)
	p.hat <- plogis(W %*% zeta.hat)

	lambda0.hat <- exp(X %*% beta0.hat)
	nu0.hat <- rep(1, n)
	p0.hat <- plogis(W %*% zeta0.hat)

	ff <- dzicmp(y, lambda.hat, nu.hat, p.hat, max = max)
	ff0 <- dzicmp(y, lambda0.hat, nu0.hat, p0.hat, max = max)

	X2 <- 2*(sum(log(ff)) - sum(log(ff0)))
	pvalue <- pchisq(X2, df = length(gamma.hat), lower.tail = FALSE)
	list(stat = X2, pvalue = pvalue)
}
