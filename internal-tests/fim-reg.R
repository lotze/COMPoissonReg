library(COMPoissonReg)

set.seed(1237)

n = 1
X = cbind(1, rnorm(n), rnorm(n))
S = cbind(1, rnorm(n))
W = cbind(1, rnorm(n))
qq = ncol(X) + ncol(S) + ncol(W)

Beta = c(1,1,1)
Gamma = c(1,2)
Zeta = c(-1,1)

d1 = ncol(X)
d2 = ncol(S)
d3 = ncol(W)

lambda = exp(X %*% Beta)
nu = exp(S %*% Gamma)
p = plogis(W %*% Zeta)

FIM = fim.zicmp.reg(X, S, W, Beta, Gamma, Zeta)

# This transformation should also yield the correct FIM for (Beta, Gamma, Zeta)
FIM.noreg = fim.zicmp(lambda, nu, p, max = 200)
J = matrix(0, 3, qq)
J[1, 1:ncol(X)] = as.numeric(exp(X %*% Beta)) * X
J[2, 1:ncol(S) + ncol(X)] = as.numeric(exp(S %*% Gamma)) * S
J[3, 1:ncol(W) + ncol(S) + ncol(X)] = as.numeric(dlogis(W %*% Zeta)) * W
t(J) %*% FIM.noreg %*% J

# Try to evaluate compute the FIM numerically as something to check against

f = function(theta, y, X, S, W) {
	Beta = theta[1:d1]
	Gamma = theta[1:d2 + d1]
	Zeta = theta[1:d3 + d1 + d2]
	sum(d.zi.compoisson(y, lambda = exp(X %*% Beta), nu = exp(S %*% Gamma), p = plogis(W %*% Zeta), log = TRUE))
}

if (FALSE) {
	R = 10000
	A = matrix(0, qq, qq)
	for (r in 1:R) {
		if (r %% 100 == 0) cat("Starting rep", r, "\n")
		y = r.zi.compoisson(n, lambda, nu, p)
		S = numDeriv::grad(f, c(Beta, Gamma, Zeta), y = y, X = X, S = S, W = W)
		A = A + S %*% t(S)
	}
	FIM.mc = A/R
}

max = 120
A = matrix(0, qq, qq)
for (y in 0:max) {
	ss = numDeriv::grad(f, c(Beta, Gamma, Zeta), y = y, X = X, S = S, W = W)
	ff = d.zi.compoisson(y, lambda, nu, p, log = FALSE)
	A = A + ss %*% t(ss) * as.numeric(ff)
}
FIM.outer = A

max = 120
A = matrix(0, qq, qq)
for (y in 0:max) {
	H = numDeriv::hessian(f, c(Beta, nu, Gamma), y = y, X = X, Z = Z)
	ff = d.zi.compoisson(y, lambda, nu, p, log = FALSE)
	# print(H)
	# print(ff)
	A = A + -H * as.numeric(ff)
}
FIM.grad = A
