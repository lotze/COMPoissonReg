library(COMPoissonReg)
library(Matrix)

x <- 1.25
X <- model.matrix(~ x)
S <- model.matrix(~ x)
W <- model.matrix(~ x)
beta <- rep(1.5, 2)
gamma <- rep(1.25, 2)
zeta <- rep(0, 2)

lambda <- exp(X %*% beta)
nu <- exp(S %*% gamma)
p <- plogis(W %*% zeta)

FIM <- fim.zicmp.reg(X, S, W = matrix(1,1,1), beta, gamma, zeta = 0, max = 100)
qq <- ncol(W) + ncol(S)
FIM[1:qq, 1:qq]

# This transformation should also yield the correct FIM for (Beta, Gamma, Zeta)
FIM.noreg <- fim.zicmp(lambda, nu, p = p, max = 200)
qq <- ncol(W) + ncol(S) + ncol(X)
J <- bdiag(
	as.numeric(exp(X %*% beta)) * X,
	as.numeric(exp(S %*% gamma)) * S,
	as.numeric(dlogis(W %*% zeta)) * W	
)

FIM.outer <- t(J) %*% FIM.noreg %*% J
FIM.outer[-c(5,6), -c(5,6)]
