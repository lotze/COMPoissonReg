library(COMPoissonReg)

set.seed(1234)

# ----- Generate the data -----
n = 400
x = runif(n, -2, 2)
X = model.matrix(~ x)
S = matrix(1, n, 1)
W = model.matrix(~ x)
beta.true = c(-0.25, 0.6)
gamma.true = -0.3
zeta.true = c(-2.5, 0.05)
lambda.true = exp(X %*% beta.true)
nu.true = exp(S %*% gamma.true)
p.true = plogis(W %*% zeta.true)
y = rzicmp(n, lambda = lambda.true, nu = nu.true, p = p.true)
dat = data.frame(y = y, x = x)

# ----- Fit ZICMP model -----
options(COMPoissonReg.optim.method = "BFGS")
zicmp.out = glm.cmp(y ~ x, formula.nu = ~ 1, formula.p = ~ x, data = dat)
print(zicmp.out)

# ----- Residuals -----
y.hat = predict(zicmp.out)
res = resid(zicmp.out, type = "quantile")
qqnorm(res); qqline(res, lty = 2, col = "red", lwd = 2)
plot(y.hat, res)

# ----- Test for equidispersion -----
equitest(zicmp.out)

# ----- Deviance -----
deviance(zicmp.out)

# ----- Bootstrap -----
boot.out = parametric.bootstrap(zicmp.out, reps = 50, report.period = 10)
plot(density(boot.out[,1])); abline(v = beta.true[1], lty = 2, lwd = 2, col = "red")
plot(density(boot.out[,2])); abline(v = beta.true[2], lty = 2, lwd = 2, col = "red")
plot(density(boot.out[,3])); abline(v = gamma.true[1], lty = 2, lwd = 2, col = "red")
plot(density(boot.out[,4])); abline(v = zeta.true[1], lty = 2, lwd = 2, col = "red")
plot(density(boot.out[,5])); abline(v = zeta.true[2], lty = 2, lwd = 2, col = "red")
