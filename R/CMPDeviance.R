CMPDeviance <- function(x,y,betahat,nuhat,leverage,max = 150){

# Computes the COM-Poisson deviances exactly

# x = matrix of size nxp where i=1,...,n and j=1,...,p
# y = col vector of size n, i=1,...,n
# betahat = COM-Poisson MLEs for beta vector
# nuhat = COM-Poisson estimate for dispersion parameter
# max = maximum number used to approximate infinite sum calculation

#### Compute optimal log likelihood value for given nu-hat value
betainit <- betahat
OptimalLogLi <- rep(0,length(y))
iterct <- rep(0,length(y))

for(i in 1:length(y)){
  # Create -logL = -logf (because considering single observation) so that we take the minimum of this function (which equals the max of logL)
    minuslogf <- function(par){-((y[i] * (x[i,] %*% par[1:length(betainit)])) - (nuhat * log(factorial(y[i]))) - log(computez(exp(x[i,] %*% par[1:length(betainit)]),nuhat,max)))}

  # Determine the MLEs
    BetaEstResult <- nlm(p=betainit,f=minuslogf)
    OptimalLogLi[i] <- -BetaEstResult$min
  }


#### Compute exact deviances


lambdahat <- exp(x %*% betahat)

OptimalLogL.mu <- (y*log(lambdahat)) - (nuhat * log(factorial(y))) - log(computez(lambdahat,nuhat,max))
OptimalLogL.y <- OptimalLogLi

d <- -2*(OptimalLogL.mu - OptimalLogL.y)

cmpdev <- d/(sqrt(1-leverage))
return(cmpdev)
}
