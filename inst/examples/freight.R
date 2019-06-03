library(COMPoissonReg)

data(freight)

options(COMPoissonReg.optim.method = "L-BFGS-B")

cmp.out = glm.cmp(broken ~ transfers, data = freight)
print(cmp.out)
coef(cmp.out)
vcov(cmp.out)
sdev(cmp.out)
nu(cmp.out)
resid(cmp.out, type = "quantile")
resid(cmp.out, type = "raw")
equitest(cmp.out)

options(COMPoissonReg.optim.method = "BFGS")

zicmp.out = glm.cmp(formula.lambda = broken ~ 1,
	formula.nu = ~ 1,
	formula.p = ~ 1,
	data = freight)
print(zicmp.out)
vcov(cmp.out)
sdev(cmp.out)
resid(zicmp.out, type = "quantile")
resid(zicmp.out, type = "raw")
equitest(zicmp.out)
