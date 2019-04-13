library(COMPoissonReg)

set.seed(1234)

# ----- CMP -----
n <- 200
x1 <- as.factor(sample(1:2, size = n, replace = TRUE))
x2 <- as.factor(sample(1:2, size = n, replace = TRUE))
X <- model.matrix(~ x1)
S <- model.matrix(~ x2)
offx <- rep(log(10), n)
offs <- rep(log(3), n)
Beta.true <- c(1, -1)
Gamma.true <- c(-1, 1)
lambda.true <- exp(X %*% Beta.true + offx)
nu.true <- exp(S %*% Gamma.true + offs)
y <- rcmp(n, lambda.true, nu.true)

out <- glm.cmp(y ~ x1 + offset(offx), formula.nu = ~ x2 + offset(offs))
print(out)

# A little simulation to make sure we're seeing the right properties
R <- 200
d <- c(length(Beta.true), length(Gamma.true))
est <- matrix(NA, R, sum(d))
for (r in 1:R) {
	cat("Rep", r, "\n")
	y <- rcmp(n, lambda.true, nu.true)
	out <- glm.cmp(y ~ x1 + offset(offx), formula.nu = ~ x2 + offset(offs))
	est[r,] <- coef(out)
}

W.stat <- numeric(R)
V.sim <- var(est)
for (r in 1:R) {
	delta <- est[r,] - c(Beta.true, Gamma.true)
	W.stat[r] <- delta %*% solve(V.sim, delta)
}

# These two should roughly match
plot(ecdf(W.stat))
curve(pchisq(x, df = sum(d)), lty = 2, lwd = 2, col = "blue", add = TRUE)

# ----- ZICMP -----
n <- 200
x1 <- as.factor(sample(1:2, size = n, replace = TRUE))
x2 <- as.factor(sample(1:2, size = n, replace = TRUE))
x3 <- as.factor(sample(1:2, size = n, replace = TRUE))
X <- model.matrix(~ x1)
S <- model.matrix(~ x2)
W <- model.matrix(~ x3)
offx <- rep(log(10), n)
offs <- rep(log(3), n)
offw <- rep(-log(3), n)
Beta.true <- c(1, -1)
Gamma.true <- c(-1, 1)
Zeta.true <- c(-1, 1)
lambda.true <- exp(X %*% Beta.true + offx)
nu.true <- exp(S %*% Gamma.true + offs)
p.true <- plogis(W %*% Gamma.true + offw)
y <- rzicmp(n, lambda.true, nu.true, p.true)

out <- glm.cmp(y ~ x1 + offset(offx),
	formula.nu = ~ x2 + offset(offs),
	formula.p = ~ x3 + offset(offw))
print(out)

# A little simulation to make sure we're seeing the right properties
R <- 200
d <- c(length(Beta.true), length(Gamma.true), length(Zeta.true))
est <- matrix(NA, R, sum(d))
for (r in 1:R) {
	cat("Rep", r, "\n")
	y <- rzicmp(n, lambda.true, nu.true, p.true)
	tryCatch({
		out <- glm.cmp(y ~ x1 + offset(offx),
			formula.nu = ~ x2 + offset(offs),
			formula.p = ~ x3 + offset(offw))
		est[r,] <- coef(out)
	}, error = function(e){
		print(e)
	})
}

W.stat <- numeric(R)
V.sim <- var(est)
for (r in 1:R) {
	delta <- est[r,] - c(Beta.true, Gamma.true, Zeta.true)
	W.stat[r] <- delta %*% solve(V.sim, delta)
}

# These two should roughly match
plot(ecdf(W.stat))
curve(pchisq(x, df = sum(d)), lty = 2, lwd = 2, col = "blue", add = TRUE)
