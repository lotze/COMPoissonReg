library(COMPoissonReg)
set.seed(1234)

# Test formulas for nu and p where only an intercept are specified. Make sure
# that the dimension is either interpreted correctly, or we throw an error.

y = rzicmp(250, lambda = 10, nu = 0.95, p = 0.05)

# ----- CMP -----
# Doesn't work
out = glm.cmp(y ~ 1)
out = glm.cmp(y ~ 1, formula.nu = ~ 1)

# Works
out = glm.cmp(y ~ 1, formula.nu = y ~ 1)

dat = data.frame(y = y)
out = glm.cmp(y ~ 1, formula.nu = ~ 1, data = dat)

# ----- ZICMP----- 
# Doesn't work
out = glm.cmp(y ~ 1, formula.nu = y ~ 1, formula.p = ~ 1)

# Works
out = glm.cmp(y ~ 1, formula.nu = y ~ 1, formula.p = y ~ 1)

dat = data.frame(y = y)
out = glm.cmp(y ~ 1, , formula.nu = ~ 1, formula.p = ~ 1, data = dat)
