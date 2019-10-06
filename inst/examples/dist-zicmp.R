library(COMPoissonReg)

# Compare ZICMP density to a histogram of draws
n = 200000
lambda = 10
nu = 0.8
p = 0.1
x = rzicmp(n, lambda, nu, p)

tt = numeric(max(x)+1)
for (i in 0:max(x)) {
	tt[i+1] = sum(x == i) / n
}
ff = dzicmp(0:max(x), lambda = lambda, nu = nu, p = p)

plot(0:max(x), tt)
points(0:max(x), ff, pch = 2, col = "blue")

FF = pzicmp(0:max(x), lambda = lambda, nu = nu, p = p)
plot(ecdf(x))
points(0:max(x), FF, pch = 2, col = "blue")

delta = min(FF)/2
qzicmp(FF - delta, lambda = lambda, nu = nu, p = p)
