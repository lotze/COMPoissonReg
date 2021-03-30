library(COMPoissonReg)

# Compare CMP density to a histogram of draws
n = 200000
lambda = 10
nu = 0.8
x = rcmp(n, lambda, nu)

tt = numeric(max(x)+1)
for (i in 0:max(x)) {
	tt[i+1] = sum(x == i) / n
}
ff = dcmp(0:max(x), lambda = lambda, nu = nu)

plot(0:max(x), tt)
points(0:max(x), ff, pch = 2, col = "blue")

# Empirical CDF vs. CDF function
FF = pcmp(0:max(x), lambda = lambda, nu = nu)
plot(ecdf(x))
points(0:max(x), FF, pch = 2, col = "blue")

delta = min(FF)/2
qcmp(FF - delta, lambda = lambda, nu = nu)

# Empirical mean vs. expected value
mean(x)
ecmp(lambda, nu)

# Empirical variance vs. variance function
var(x)
vcmp(lambda, nu)
