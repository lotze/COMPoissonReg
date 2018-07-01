library(COMPoissonReg)

# In the examples below, distn of y has nothing to do with x, so correct
# regression coefficients are zero. But the large counts might screw us
# up numerically, if the means are not accounting for them.

options(COMPoissonReg.optim.method = "BFGS")

y <- c(320, 600, 312, 549, 246, 65, 41, 52)
x <- c(0, 0, 1, 1, 1, 1, 1, 1)
glm.cmp(y ~ x, data=data.frame(x,y))

y <- c(167, 566, 64, 98, 33, 47, 40, 27)
x1 <- c(0, 0, 1, 1, 0, 0, 1, 1)
x2 <- c(0, 0, 0, 0, 1, 1, 1, 1)
glm.cmp(y ~ x1 + x2, data=data.frame(y,x1,x2))
glm.cmp(y ~ x1 + x2, data=data.frame(y,x1,x2), beta.init = c(4,0,0), gamma.init = log(0.25))
glm.cmp(y ~ 1, data=data.frame(y,x1,x2), beta.init = 6, gamma.init = log(0.25))

y <- rcmp(length(x1), lambda = exp(2.25), nu = 0.4)
glm.cmp(y ~ 1, data=data.frame(y,x1,x2), beta.init = 6, gamma.init = log(1))

x <- sample(c(0,1), size = 800, replace = TRUE)
y <- rcmp(length(x), lambda = exp(-0.25 + 6*x), nu = 0.95)
glm.cmp(y ~ x, data=data.frame(x,y))

# Very small nu causes numerical issues with some functions
tryCatch({
	# Exact and approx z-function evaluate to Inf
	z_exact(lambda = 10, nu = 0.001, take_log = TRUE)
	z_approx(lambda = 10, nu = 0.001, take_log = TRUE)

	# Draws are always the same (somewhat large) number. This occurs because the
	# adaptive z-function evaluates to infinity before it reaches any significant
	# mass in the CMP distribution.
	rcmp(20, lambda = 10, nu = 0.001)

	# Expected value becomes excessively large as nu -> 0, and our function
	# evaluates to NaN
	cmp_expected_value(lambda = 10, nu = 0.001)
}, error = function(e) {
	print(e)
})