InfoMatrix <- function(I.beta, I.betanu, I.nu)
{
	# Put the components/blocks that comprise the full information matrix together
	I1 <- cbind(I.beta,I.betanu)
	I2 <- cbind(t(I.betanu),I.nu)

	# Create the information matrix
	I <- rbind(I1,I2)
	return(I)
}

InfoMatrix.beta <- function(x,weight)
{
	#create vector of ones
	if(is.matrix(x)==TRUE || is.data.frame(x) == TRUE) {onevec <- rep(1,length(x[,1]))} else onevec <- rep(1,length(x)) 

	#create real X matrix, namely where 1st col is vector of 1s to incorporate beta0 effect
	newx <- cbind(onevec,x) 
	xmat <- as.matrix(newx)

	# Create the I(beta) block to the information matrix
	W <- diag(weight)
	I <- t(xmat) %*% W %*% xmat
	return(I)
}

InfoMatrix.betanu <- function(x,beta,nu,max)
{
	# create vector of ones
	if(is.matrix(x)==TRUE || is.data.frame(x) == TRUE) {onevec <- rep(1,length(x[,1]))} else onevec <- rep(1,length(x)) 

	# create real X matrix, namely where 1st col is vector of 1s to incorporate beta0 effect
	newx <- cbind(onevec,x) 
	xmat <- as.matrix(newx)

	# Compute the components needed to determine the beta,nu block of the information matrix
	ans1 <- computez.prodjlogj(exp(xmat %*% beta), nu, max)
	ans2 <- computez.prodlogj(exp(xmat %*% beta), nu, max)
	ans3 <- computez.prodj(exp(xmat %*% beta), nu, max)
	ans4 <- computez(exp(xmat %*% beta), nu, max)

	ans <- (ans1/ans4) - ((ans2/ans4) * (ans3/ans4))

	# Compute the beta,nu block of the information matrix
	I <- t(xmat) %*% ans
	return(I)
}

InfoMatrix.nu <- function(x,beta,nu,max)
{
	#create vector of ones
	if(is.matrix(x)==TRUE || is.data.frame(x) == TRUE) {onevec <- rep(1,length(x[,1]))} else onevec <- rep(1,length(x)) 

	#create real X matrix, namely where 1st col is vector of 1s to incorporate beta0 effect
	newx <- cbind(onevec,x) 
	xmat <- as.matrix(newx)

	# Compute the components needed to determine the information component due to nu
	ans1 <- computez.prodlogj2(exp(xmat %*% beta), nu, max)
	ans2 <- computez.prodlogj(exp(xmat %*% beta), nu, max)
	ans3 <- computez(exp(xmat %*% beta), nu, max)

	# Compute the information associated with nu
	ans <- ans1/ans3 - ((ans2/ans3)^2)
	return(sum(ans))
}
