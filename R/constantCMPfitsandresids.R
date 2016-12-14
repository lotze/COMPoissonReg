constantCMPfitsandresids <- function(betahat,nuhat,x,y=0){

#Determine estimated lambdahat, fit, and residuals
   lambdahat <- exp(x %*% betahat)
   fit <- lambdahat^(1/nuhat) - ((nuhat -1)/(2*nuhat))

   resid <- y - fit

return(list(fit=fit,resid=resid))
}

