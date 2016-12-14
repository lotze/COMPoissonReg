ComputeBetasAndNuHat <- function(x, y, betainit, nuinit, max)
{
	# Uses nlminb to solve for the MLE estimates for betas and nu
	
	# x = matrix of size nxp where i=1,...,n and j=1,...,p
	# y = col vector of size n, i=1,...,n
	# betainit = initial vector of betas, b0_1, ..., b0_p
	# nuinit = initial nu value

	# Create -logL so that we take the minimum of this function (which equals the max of logL)
	minusloglike <- function(par) {
		beta <- par[1:length(betainit)]
		nu <- par[length(betainit)+1]
		z <- computez(exp(x %*% beta), nu, max)
   		-sum((y * (x %*% beta)) - (nu * lgamma(y+1)) - log(z))
   	}

	# Determine the MLEs
	BetaNuEst <- nlminb(start=c(betainit,nuinit),minusloglike,lower = c(rep(-Inf,length(betainit)),0), upper = c(rep(Inf,length(betainit)),Inf))

	if (BetaNuEst$convergence != 0) {
		# Determine the MLEs
		BetaNuEst <- optim(par=c(betainit,nuinit),fn=minusloglike,control=list(maxit=1000))#$par
	}

	return(BetaNuEst)
}

