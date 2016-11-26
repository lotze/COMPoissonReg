library(COMPoissonReg)

set.seed(1235)

# ----- Generate the data -----
n <- 400
x <- runif(n, 1, 4)
X <- model.matrix(~ x)
beta.true <- c(1, 0.5)
zeta.true <- c(0.25, -1)
lambda.true <- exp(X %*% beta.true)
nu.true <- 0.75
p.true <- plogis(X %*% zeta.true)

y <- numeric(n)
for (i in 1:n) {
	y[i] <- r.zi.compoisson(1, lambda = lambda.true[i], nu = nu.true, p = p.true[i])
}
dat <- data.frame(y = y, x = x)

# ----- Fit ZICMP model -----
zicmp.out <- zicmp(y ~ x, formula.nu = ~ 1, formula.p = ~ x, data = dat)
print(zicmp.out)

# ----- Residuals -----
y.hat <- predict(zicmp.out)
res <- resid(zicmp.out, type = "quantile")
qqnorm(res); qqline(res, lty = 2, col = "red", lwd = 2)
plot(y.hat, res)
