LRT <- function(x,y,betahat0,betahat,nuhat,max){
# This function computes the -2logLRT value and associated p-value for significance

#create vector of ones
  if(is.matrix(x)==TRUE || is.data.frame(x) == TRUE) {onevec <- rep(1,length(x[,1]))} else onevec <- rep(1,length(x)) 

#create real X matrix, namely where 1st col is vector of 1s to incorporate beta0 effect
  newx <- cbind(onevec,x) 
  xmat <- as.matrix(newx)

# Compute the test statistic
  teststat <- -2*((t(y)%*%(xmat %*% betahat0)) - sum(log(factorial(y))) - sum(exp(xmat %*% betahat0)) - 
              ((t(y)%*%(xmat %*% betahat)) - (nuhat*sum(log(factorial(y)))) - 
              sum(log(computez(exp(xmat %*% betahat),nuhat,max)))))

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

	for (i in 1:n) {
		ff[i] <- d.zi.compoisson(y[i], lambda.hat[i], nu.hat[i], p.hat[i], max)
		ff0[i] <- d.zi.compoisson(y[i], lambda0.hat[i], nu0.hat[i], p0.hat[i], max)
	}

	X2 <- 2*(sum(log(ff)) - sum(log(ff0)))
	pvalue <- pchisq(X2, df = 1, lower.tail = FALSE)
	list(stat = X2, pvalue = pvalue)
}
