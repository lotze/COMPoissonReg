library(COMPoissonReg)

data(freight)

cmp.out <- glm.cmp(formula = broken ~ transfers, data = freight)
print(cmp.out)
coef(cmp.out)
nu(cmp.out)
resid(cmp.out, type = "quantile")
resid(cmp.out, type = "raw")

zicmp.out <- glm.zicmp(formula.lambda = broken ~ transfers, 
#	formula.nu = ~ transfers,
	formula.p = ~ transfers,
	data = freight)
print(zicmp.out)

resid(zicmp.out)
resid(zicmp.out, type = "quantile")

