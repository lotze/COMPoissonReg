CMPStdErrors <-
function(x,beta,nu,max=100){

   weight <- weights(x,beta,nu,max)
   Ibeta <- InfoMatrix.beta(x,weight)
   Ibetanu <- InfoMatrix.betanu(x,beta,nu,max)
   Inu <- InfoMatrix.nu(x,beta,nu,max)
   info <- InfoMatrix(Ibeta,Ibetanu,Inu)

   # If we want to put a regression on nu, we can get Std Errors like this
   # via ZICMP FIM. [Do we want to support that regression?]
   # qq <- ncol(X) + ncol(S)
   # FF <- fim.zicmp.reg(X, S, W = matrix(1,1,1), beta, gamma, zeta = 0, max = max)
   # info <- FF[1:qq, 1:qq]

   CMP.SEs <- sqrt(diag(solve(info)))

return(CMP.SEs)

}

