library(COMPoissonReg)

set.seed(1234)

# ----- Generate the data -----
n <- 400
x <- runif(n, 1, 4)
X <- model.matrix(~ x)
beta.true <- c(1, 0.5)
lambda.true <- exp(X %*% beta.true)
nu.true <- 0.75

y <- numeric(n)
for (i in 1:n) {
	y[i] <- rcom(1, lambda = lambda.true[i], nu = nu.true)
}
dat <- data.frame(y = y, x = x)

# ----- Fit ZICMP model -----
cmp.out <- cmp(y ~ x, data = dat)
print(cmp.out)

# ----- Residuals -----
res <- resid(cmp.out, type = "quantile")
y.hat <- predict(cmp.out)
qqnorm(res); qqline(res, lty = 2, col = "red", lwd = 2)
plot(y.hat, res)

# ----- Test for equidispersion -----
chisq(cmp.out)
pval(cmp.out)

# ----- Deviance -----
deviance(cmp.out)

# ----- Bootstrap -----
boot.out <- parametric_bootstrap(cmp.out, reps = 50, report.period = 10)
hist(boot.out[,1])
hist(boot.out[,2])
hist(boot.out[,3])
