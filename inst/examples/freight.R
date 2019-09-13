library(COMPoissonReg)

data(freight)

options(COMPoissonReg.optim.method = "L-BFGS-B")

cmp_out = glm_cmp(broken ~ transfers, data = freight)
print(cmp_out)
coef(cmp_out)
vcov(cmp_out)
sdev(cmp_out)
nu(cmp_out)
resid(cmp_out, type = "quantile")
resid(cmp_out, type = "raw")
equitest(cmp_out)

options(COMPoissonReg.optim.method = "BFGS")

zicmp_out = glm_cmp(formula_lambda = broken ~ 1,
	formula_nu = ~ 1,
	formula_p = ~ 1,
	data = freight)
print(zicmp_out)
vcov(cmp_out)
sdev(cmp_out)
resid(zicmp_out, type = "quantile")
resid(zicmp_out, type = "raw")
equitest(zicmp_out)
