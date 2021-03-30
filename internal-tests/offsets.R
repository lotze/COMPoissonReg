library(COMPoissonReg)
set.seed(1234)

# ----- CMP -----
n = 200
x1 = as.factor(sample(1:2, size = n, replace = TRUE))
x2 = as.factor(sample(1:2, size = n, replace = TRUE))
X = model.matrix(~ x1)
S = model.matrix(~ x2)
offx = rep(log(10), n)
offs = rep(log(3), n)
beta.true = c(1, -1)
gamma.true = c(-1, 1)
lambda.true = exp(X %*% beta.true + offx)
nu.true = exp(S %*% gamma.true + offs)
y = rcmp(n, lambda.true, nu.true)

dat = data.frame(y, x1, x2, offx, offs)
out = glm.cmp(y ~ x1 + offset(offx), formula.nu = ~ x2 + offset(offs), data = dat)
print(out)
predict(out)

# Predict with new data, and compare to manual computation
n_new = 50
dat.new = data.frame(
	x1 = as.factor(sample(1:2, size = n_new, replace = TRUE)),
	x2 = as.factor(sample(1:2, size = n_new, replace = TRUE)),
	offx = rep(log(10), n_new),
	offs = rep(log(3), n_new))
predict(out, newdata = dat.new)

X.new = model.matrix(~ x1 + offset(offx), data = dat.new)
S.new = model.matrix(~ x2 + offset(offs), data = dat.new)
lambda.new = exp(X.new %*% out$beta + dat.new$offx)
nu.new = exp(S.new %*% out$gamma + dat.new$offs)
ecmp(lambda.new, nu.new)

# A little simulation to check the properties
R = 200
d = c(length(beta.true), length(gamma.true))
est = matrix(NA, R, sum(d))
for (r in 1:R) {
	cat("Rep", r, "\n")
	y = rcmp(n, lambda.true, nu.true)
	dat.sim = data.frame(y, x1, x2, offx, offs)
	out = glm.cmp(y ~ x1 + offset(offx), formula.nu = ~ x2 + offset(offs),
		data = dat.sim)
	est[r,] = coef(out)
}

W.stat = numeric(R)
V.sim = var(est)
for (r in 1:R) {
	delta = est[r,] - c(beta.true, gamma.true)
	W.stat[r] = delta %*% solve(V.sim, delta)
}

# These two should roughly match
plot(ecdf(W.stat))
curve(pchisq(x, df = sum(d)), lty = 2, lwd = 2, col = "blue", add = TRUE)

# ----- ZICMP -----
n = 200
x1 = as.factor(sample(1:2, size = n, replace = TRUE))
x2 = as.factor(sample(1:2, size = n, replace = TRUE))
x3 = as.factor(sample(1:2, size = n, replace = TRUE))
X = model.matrix(~ x1)
S = model.matrix(~ x2)
W = model.matrix(~ x3)
offx = rep(log(10), n)
offs = rep(log(3), n)
offw = rep(-log(3), n)
beta.true = c(1, -1)
gamma.true = c(-1, 1)
zeta.true = c(-1, 1)
lambda.true = exp(X %*% beta.true + offx)
nu.true = exp(S %*% gamma.true + offs)
p.true = plogis(W %*% zeta.true + offw)
y = rzicmp(n, lambda.true, nu.true, p.true)

dat = data.frame(y, x1, x2, x3, offx, offs, offw)
out = glm.cmp(y ~ x1 + offset(offx),
	formula.nu = ~ x2 + offset(offs),
	formula.p = ~ x3 + offset(offw), data = dat)
print(out)

# Predict with new data, and compare to manual computation
n_new = 50
dat.new = data.frame(
	x1 = as.factor(sample(1:2, size = n_new, replace = TRUE)),
	x2 = as.factor(sample(1:2, size = n_new, replace = TRUE)),
	x3 = as.factor(sample(1:2, size = n_new, replace = TRUE)),
	offx = rep(log(10), n_new),
	offs = rep(log(3), n_new),
	offw = rep(log(2), n_new))
predict(out, newdata = dat.new)

X.new = model.matrix(~ x1 + offset(offx), data = dat.new)
S.new = model.matrix(~ x2 + offset(offs), data = dat.new)
W.new = model.matrix(~ x3 + offset(offw), data = dat.new)
lambda.new = exp(X.new %*% out$beta + dat.new$offx)
nu.new = exp(S.new %*% out$gamma + dat.new$offs)
p.new = plogis(W.new %*% out$zeta + dat.new$offw)
ezicmp(lambda.new, nu.new, p.new)

# A little simulation to check the properties
R = 200
d = c(length(beta.true), length(gamma.true), length(zeta.true))
est = matrix(NA, R, sum(d))
for (r in 1:R) {
	cat("Rep", r, "\n")
	y = rzicmp(n, lambda.true, nu.true, p.true)
	dat.sim = data.frame(y, x1, x2, x3, offx, offs, offw)
	tryCatch({
		out = glm.cmp(y ~ x1 + offset(offx),
			formula.nu = ~ x2 + offset(offs),
			formula.p = ~ x3 + offset(offw), data = dat.sim)
		est[r,] = coef(out)
	}, error = function(e){
		print(e)
	})
}

W.stat = numeric(R)
V.sim = var(est)
for (r in 1:R) {
	delta = est[r,] - c(beta.true, gamma.true, zeta.true)
	W.stat[r] = delta %*% solve(V.sim, delta)
}

# These two should roughly match
plot(ecdf(W.stat))
curve(pchisq(x, df = sum(d)), lty = 2, lwd = 2, col = "blue", add = TRUE)
