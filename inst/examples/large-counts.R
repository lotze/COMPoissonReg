library(COMPoissonReg)
set.seed(1234)

n <- 400
X <- matrix(1, n, 1)
beta.true <- log(2)
Xbeta.true <- X %*% beta.true
lambda.true <- exp(Xbeta.true)
gamma.true <- log(0.2)
nu.true <- exp(gamma.true)
y <- rcmp(n, lambda.true, nu.true)
dat <- data.frame(x = X[,1], y = y)

options(COMPoissonReg.optim.method = "L-BFGS-B")
glm.out <- glm.cmp(formula.lambda = y ~ 1,
	formula.nu = ~1,
	data = dat,
	beta.init = beta.true, gamma.init = log(nu.true))
print(glm.out)

res <- residuals(glm.out, type = "quantile")
qqnorm(res); abline(c(0,1), lty = 2, col = "red")
