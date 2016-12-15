library(COMPoissonReg)

data(couple)

zicmp.out <- glm.zicmp(formula.lambda = UPB ~ EDUCATION + ANXIETY,
	formula.nu = ~ 1,
	formula.p = ~ EDUCATION + ANXIETY,
	data = couple)
print(zicmp.out)

equitest(zicmp.out)

res <- resid(zicmp.out, type = "quantile")
qqnorm(res); qqline(res, col = "red", lwd = 2, lty = 2)
