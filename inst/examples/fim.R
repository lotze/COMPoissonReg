library(COMPoissonReg)

lambda <- exp(10.57)
nu <- 4.45
p <- 0.001

FIM1 <- fim.cmp(lambda, nu)
FIM2 <- fim.zicmp(lambda, nu, p = 0.01)
round(FIM1, 6)
round(FIM2, 6)

FIM1 <- fim.zicmp(lambda, nu, p)
FIM2 <- fim.zicmp.old(lambda, nu, p, max = 10000)
round(FIM1, 6)
round(FIM2, 6)

n <- 10
X <- matrix(1, n, 1)
S <- matrix(1, n, 1)
W <- matrix(1, n, 1)
beta <- log(lambda)
gamma <- log(nu)
zeta <- qlogis(p)

FIM1 <- fim.zicmp.reg(X, S, W, beta, gamma, zeta)
FIM2 <- fim.zicmp.reg.old(X, S, W, beta, gamma, zeta, max = 10000)
round(FIM1, 6)
round(FIM2, 6)

