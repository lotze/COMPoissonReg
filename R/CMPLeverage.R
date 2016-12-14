CMPLeverage <- function(x,y,betahat,nuhat,max){

# 1) to code the W matrix  (diagonal matrix with Var(Y_i) )

   W <- diag(weights(x,betahat,nuhat,max))

#    and X matrix (in Appendix)

   E.y <- computez.prodj(exp(x %*% betahat),nuhat,max)/computez(exp(x %*% betahat),nuhat,max)
   E.logfacty <- computez.prodlogj(exp(x %*% betahat),nuhat,max)/computez(exp(x %*% betahat),nuhat,max)
   extravec <- (-log(factorial(y)) + E.logfacty)/(y - E.y)
   curlyX.mat <- cbind(x,extravec)

# 2) to compute H using eq (12)  on p. 11

   H1 <- t(curlyX.mat) %*% sqrt(W)
   H2 <- solve(t(curlyX.mat) %*% W %*% curlyX.mat)
   H <- t(H1) %*% H2 %*% H1
   diagH <- diag(H)
return(diagH)
}

