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

nu.hat <- cmp.out$nu
beta.hat <- cmp.out$coefficients
beta0.hat <- cmp.out$glm_coefficients
LRT(X[,-1], y, beta0.hat, beta.hat, nu.hat, max = 1000)

z <- computez(lambda.hat, nu.hat, max = 100)
z0 <- computez(lambda0.hat, nu = 1, max = 100)
lambda.hat <- exp(X %*% beta.hat)
lambda0.hat <- exp(X %*% beta0.hat)

ff <- numeric(n)
ff0 <- numeric(n)
for (i in 1:n) {
	ff[i] <- dcom(y[i], lambda.hat[i], nu.hat, z = z[i])
	ff0[i] <- dcom(y[i], lambda0.hat[i], 1, z = z0[i])
}
X2 <- 2*(sum(log(ff)) - sum(log(ff0)))
pvalue <- pchisq(X2, df = 1, lower.tail = FALSE)


	