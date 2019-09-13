library(COMPoissonReg)

data(couple)

cmp_out = glm_cmp(formula_lambda = UPB ~ EDUCATION + ANXIETY,
	formula_nu = ~ 1,
	formula_p = ~ EDUCATION + ANXIETY,
	data = couple)
print(cmp_out)

vcov(cmp_out)
sdev(cmp_out)

equitest(cmp_out)

res = resid(cmp_out, type = "quantile")
qqnorm(res); qqline(res, col = "red", lwd = 2, lty = 2)
