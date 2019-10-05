library(COMPoissonReg)
set.seed(1234)

n = 400
X = matrix(1, n, 1)
beta.true = log(2)
Xbeta.true = X %*% beta.true
lambda.true = exp(Xbeta.true)
gamma.true = log(0.2)
nu.true = exp(gamma.true)
y = rcmp(n, lambda.true, nu.true)
dat = data.frame(x = X[,1], y = y)

options(COMPoissonReg.optim.method = "L-BFGS-B")
glm.out = glm.cmp(formula.lambda = y ~ 1,
	formula.nu = ~1,
	data = dat,
	beta.init = beta.true, gamma.init = log(nu.true))
print(glm.out)

res = residuals(glm.out, type = "quantile")
qqnorm(res); abline(c(0,1), lty = 2, col = "red")


# Some combinations of parameters move the mass of CMP to very large values.
# This can lead to functions like rcmp taking a huge amount of computation
# and memory if we're not careful. We set a global ymax parameter to protect
# against this.

n = 5
offx = rep(5, n)
offs = rep(-5, n)
lambda.true = rep(3, n)
nu.true = rep(0.1, n)

options(COMPoissonReg.ymax = 1e3)
y = rcmp(n, lambda.true, nu.true)

options(COMPoissonReg.ymax = 1e4)
y = rcmp(n, lambda.true, nu.true)

options(COMPoissonReg.ymax = 1e5)
y = rcmp(n, lambda.true, nu.true)

options(COMPoissonReg.ymax = 1e6)
y = rcmp(n, lambda.true, nu.true)

# Make sure the offset term is interpreted correctly
n = 256
offx = rep(-1, n)
offw = rep(-1, n)
X = matrix(1, n, 1)
beta.true = log(2)
Xbeta.true = X %*% beta.true
lambda.true = exp(Xbeta.true + offx)
gamma.true = log(1)
nu.true = exp(gamma.true + offw)
y = rcmp(n, lambda.true, nu.true)
dat = data.frame(x = 1, y = y, offx = offx, offw = offw)

glm.out = glm.cmp(
	formula.lambda = y ~ 1 + offset(offx),
	formula.nu = ~1 + offset(offw),
	data = dat,
	beta.init = beta.true, gamma.init = gamma.true)
print(glm.out)

# Run a small simulation to check the distribution of beta and gamma estimators
R = 200
coef.sim = matrix(NA, R, 2)
for (r in 1:R) {
	if (r %% 25 == 0) { cat("Starting rep", r, "\n") }
	y = rcmp(n, lambda.true, nu.true)
	dat = data.frame(x = 1, y = y, offx = offx)
	glm.out = glm.cmp(
		formula.lambda = y ~ 1 + offset(offx),
		formula.nu = ~1 + offset(offw),
		data = dat,
		beta.init = beta.true, gamma.init = gamma.true)
	coef.sim[r,] = coef(glm.out)
}
V.sim = var(coef.sim)

w = numeric(R)
for (r in 1:R) {
	delta = coef.sim[r,] - c(beta.true, gamma.true)
	w[r] = t(delta) %*% solve(V.sim, delta)
}
plot(ecdf(w))
curve(pchisq(x, df = 2), lwd = 2, col = "blue", add = TRUE)
