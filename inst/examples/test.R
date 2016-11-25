library(COMPoissonReg)

data(freight)

cmp.out <- cmp(formula = broken ~ transfers, data = freight)
print(cmp.out)
coef(cmp.out)
nu(cmp.out)

zicmp.out <- zicmp(formula.lambda = broken ~ transfers, 
	formula.nu = ~ transfers,
	formula.p = ~ transfers,
	data = freight)
print(zicmp.out)

zicmp.out$opt.res$message
