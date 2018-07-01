# In this case, distn of y has nothing to do with x, so correct regression coefficients are zero
# But the large counts might screw us up numerically, if the means are not accounting for them

# Another problem is with rcmp: Our usual truncation of the sum at 100 cannot handle generation
# of large counts. Also, the R code is currently way to slow to generate large counts with a 
# large max.

y <- c(320, 600, 312, 549, 246, 65, 41, 52)
x <- c(0, 0, 1, 1, 1, 1, 1, 1)
glm.cmp(y ~ x, data=data.frame(x,y))

y <- c(167, 566, 64, 98, 33, 47, 40, 27)
x1 <- c(0, 0, 1, 1, 0, 0, 1, 1)
x2 <- c(0, 0, 0, 0, 1, 1, 1, 1)
glm.cmp(y ~ x1 + x2, data=data.frame(y,x1,x2))
glm.cmp(y ~ x1 + x2, data=data.frame(y,x1,x2), beta.init = c(4,0,0), gamma.init = log(0.25))
glm.cmp(y ~ 1, data=data.frame(y,x1,x2), beta.init = 6, gamma.init = log(0.25))

options(COMPoissonReg.optim.method = "BFGS")
options(COMPoissonReg.optim.method = "L-BFGS-B")
options(COMPoissonReg.optim.method = "Nelder-Mead")

tryCatch({
	y <- rcmp(1, lambda = exp(5.25), nu = 0.4)
}, error = function(e) {
	print(e)
})

y <- rcmp(length(x1), lambda = exp(2.25), nu = 0.4)
glm.cmp(y ~ 1, data=data.frame(y,x1,x2), beta.init = 6, gamma.init = log(1))

x <- sample(c(0,1), size = 800, replace = TRUE)
y <- rcmp(length(x), lambda = exp(-0.25 + 6*x), nu = 0.95)
glm.cmp(y ~ x, data=data.frame(x,y))
