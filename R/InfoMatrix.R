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
	# Create the I(beta) block to the information matrix
	W <- diag(weight)
	I <- t(x) %*% W %*% x
	return(I)
}

InfoMatrix.betanu <- function(x,beta,nu,max)
{
	# Compute the components needed to determine the beta,nu block of the information matrix
	ans1 <- computez.prodjlogj(exp(x %*% beta), nu, max)
	ans2 <- computez.prodlogj(exp(x %*% beta), nu, max)
	ans3 <- computez.prodj(exp(x %*% beta), nu, max)
	ans4 <- computez(exp(x %*% beta), nu, max)

	ans <- (ans1/ans4) - ((ans2/ans4) * (ans3/ans4))

	# Compute the beta,nu block of the information matrix
	I <- t(x) %*% ans
	return(I)
}

InfoMatrix.nu <- function(x,beta,nu,max)
{
	# Compute the components needed to determine the information component due to nu
	ans1 <- computez.prodlogj2(exp(x %*% beta), nu, max)
	ans2 <- computez.prodlogj(exp(x %*% beta), nu, max)
	ans3 <- computez(exp(x %*% beta), nu, max)

	# Compute the information associated with nu
	ans <- ans1/ans3 - ((ans2/ans3)^2)
	return(sum(ans))
}
