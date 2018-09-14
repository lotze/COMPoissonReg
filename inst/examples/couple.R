library(COMPoissonReg)

data(couple)

cmp.out <- glm.cmp(formula.lambda = UPB ~ EDUCATION + ANXIETY,
	formula.nu = ~ 1,
	formula.p = ~ EDUCATION + ANXIETY,
	data = couple)
print(cmp.out)

round(vcov(cmp.out, use.fim = FALSE), 6)
round(vcov(cmp.out, use.fim = TRUE), 6)
round(vcov(cmp.out, use.fim = TRUE, mc.reps = 1000), 6)

sdev(cmp.out, use.fim = FALSE)
sdev(cmp.out, use.fim = TRUE)
sdev(cmp.out, use.fim = TRUE, mc.reps = 1000)

equitest(cmp.out)

res <- resid(cmp.out, type = "quantile")
qqnorm(res); qqline(res, col = "red", lwd = 2, lty = 2)
