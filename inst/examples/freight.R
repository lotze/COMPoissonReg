library(COMPoissonReg)

data(freight)

cmp.out <- glm.cmp(broken ~ transfers, data = freight)
print(cmp.out)
coef(cmp.out)
vcov(cmp.out, fim = FALSE)
vcov(cmp.out, fim = TRUE)
nu(cmp.out)
resid(cmp.out, type = "quantile")
resid(cmp.out, type = "raw")
equitest(cmp.out)

options(COMPoissonReg.optim.method = "BFGS")

zicmp.out <- glm.cmp(formula.lambda = broken ~ 1,
	formula.nu = ~ 1,
	formula.p = ~ 1,
	data = freight)
print(zicmp.out)
vcov(cmp.out, fim = FALSE)
vcov(cmp.out, fim = TRUE)
resid(zicmp.out)
resid(zicmp.out, type = "quantile")
equitest(zicmp.out)
