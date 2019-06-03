library(COMPoissonReg)

# A small simulation
n = 150
X = cbind(int = 1, x = rnorm(n, 0, 0.5))
S = cbind(int = 1, s = rnorm(n, 0, 0.5))
W = cbind(int = 1, w = rnorm(n, 0, 0.5))
offx = rep(0.5, n)
offs = rep(-0.1, n)
offw = rep(-2, n)
beta.true = c(-0.75, 1.5)
gamma.true = c(0, -1)
zeta.true = c(-1, 0.25)
lambda.true = exp(X %*% beta.true + offx)
nu.true = exp(S %*% gamma.true + offs)
p.true = plogis(W %*% zeta.true + offw)

R = 200
res = matrix(NA, R, 6)
for (r in 1:R) {
	if (r %% 1 == 0) { cat("Starting rep", r, "\n") }
	y = rzicmp(n, lambda.true, nu.true, p.true)
	dat = data.frame(y = y, x = X[,2], s = S[,2], w = W[,2],
		offx = offx, offs = offs, offw = offw)
	tryCatch({
		out = glm.cmp(
			formula.lambda = y ~ x + offset(offx),
			formula.nu = ~ s + offset(offs),
			formula.p = ~ w + offset(offw),
			data = dat)
		res[r,] = coef(out)
	}, error = function(e) {
		print(e)
	})
}

ww = numeric(R)
V.mc = var(res, na.rm = TRUE)
for (r in 1:R) {
	delta = res[r,] - c(beta.true, gamma.true, zeta.true)
	ww[r] = t(delta) %*% solve(V.mc, delta)
}

plot(ecdf(ww))
curve(pchisq(x, df = ncol(res)), add = TRUE, lty = 2, lwd = 2, col = "red")

plot(density(ww, na.rm = TRUE))
curve(dchisq(x, df = ncol(res)), add = TRUE, lty = 2, lwd = 2, col = "red")
