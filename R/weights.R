weights <- function(x,beta,nu,max){

# Compute the parts that comprise the weight functions
  w1 <- computez.prodj2(exp(x %*% beta),nu,max)
  w2 <- computez.prodj(exp(x %*% beta),nu,max)
  w3 <- computez(exp(x %*% beta),nu,max)

  Ey2 <- w1/w3
  E2y <- (w2/w3)^2

# Determine the weights
  weight <- Ey2 - E2y

return(weight)
}

