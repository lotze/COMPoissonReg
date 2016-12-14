library(COMPoissonReg)

set.seed(1235)

# ----- Generate the data -----
n <- 400
x <- runif(n, 1, 4)
X <- model.matrix(~ x)
S <- matrix(1, n, 1)
W <- model.matrix(~ x)
beta.true <- c(1, 2)
gamma.true <- 1
zeta.true <- c(0.05, -1)
lambda.true <- exp(X %*% beta.true)
nu.true <- exp(S %*% gamma.true)
p.true <- plogis(W %*% zeta.true)
y <- rzicmp(n, lambda = lambda.true, nu = nu.true, p = p.true)
dat <- data.frame(y = y, x = x)

# ----- Fit ZICMP model -----
zicmp.out <- glm.zicmp(y ~ x, formula.nu = ~ 1, formula.p = ~ x, data = dat)
print(zicmp.out)

# ----- Residuals -----
y.hat <- predict(zicmp.out)
res <- resid(zicmp.out, type = "quantile")
qqnorm(res); qqline(res, lty = 2, col = "red", lwd = 2)
plot(y.hat, res)

# ----- Test for equidispersion -----
chisq(zicmp.out)
pval(zicmp.out)

# ----- Deviance -----
deviance(zicmp.out)

# ----- Bootstrap -----
boot.out <- parametric_bootstrap(zicmp.out, reps = 50, report.period = 10)
hist(boot.out[,1])
hist(boot.out[,2])
hist(boot.out[,3])
hist(boot.out[,4])
hist(boot.out[,5])
