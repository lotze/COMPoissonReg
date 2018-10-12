library(COMPoissonReg)

data(couple)

cmp.out <- glm.cmp(formula.lambda = UPB ~ EDUCATION + ANXIETY,
	formula.nu = ~ 1,
	formula.p = ~ EDUCATION + ANXIETY,
	data = couple)
print(cmp.out)

vcov(cmp.out)
sdev(cmp.out)

equitest(cmp.out)

res <- resid(cmp.out, type = "quantile")
qqnorm(res); qqline(res, col = "red", lwd = 2, lty = 2)
