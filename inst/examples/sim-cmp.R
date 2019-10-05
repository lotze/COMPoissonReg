library(COMPoissonReg)

set.seed(1234)

# ----- Generate the data -----
n = 400
x = runif(n, 1, 4)
X = model.matrix(~ x)
beta.true = c(1, 0.5)
lambda.true = exp(X %*% beta.true)
nu.true = 0.75
y = rcmp(n, lambda = lambda.true, nu = nu.true)
dat = data.frame(y = y, x = x)

# ----- Fit ZICMP model -----
cmp.out = glm.cmp(y ~ x, data = dat)
print(cmp.out)

# ----- Residuals -----
res = resid(cmp.out, type = "quantile")
y.hat = predict(cmp.out)
qqnorm(res); qqline(res, lty = 2, col = "red", lwd = 2)
plot(y.hat, res)

# ----- Test for equidispersion -----
equitest(cmp.out)

# ----- Deviance -----
deviance(cmp.out)

# ----- Leverage -----
lev = leverage(cmp.out)
plot(y.hat, lev)

# ----- Bootstrap -----
boot.out = parametric.bootstrap(cmp.out, reps = 50, report.period = 10)
plot(density(boot.out[,1])); abline(v = beta.true[1], lty = 2, lwd = 2, col = "red")
plot(density(boot.out[,2])); abline(v = beta.true[2], lty = 2, lwd = 2, col = "red")
plot(density(boot.out[,3])); abline(v = gamma.true[1], lty = 2, lwd = 2, col = "red")
