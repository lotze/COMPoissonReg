library(COMPoissonReg)

data(freight)

cmp.out <- glm.cmp(formula = broken ~ transfers, data = freight)
print(cmp.out)
coef(cmp.out)
nu(cmp.out)
resid(cmp.out, type = "quantile")
resid(cmp.out, type = "raw")
equitest(cmp.out)

zicmp.out <- glm.zicmp(formula.lambda = broken ~ 1,
	formula.nu = ~ 1,
	formula.p = ~ 1,
	data = freight)
print(zicmp.out)

resid(zicmp.out)
resid(zicmp.out, type = "quantile")
equitest(zicmp.out)
