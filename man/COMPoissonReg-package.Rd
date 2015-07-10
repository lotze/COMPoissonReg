\name{COMPoissonReg-package}

\alias{COM-PoissonReg-package}
\alias{COM-PoissonReg}
\alias{COMPoissonReg}

\docType{package}

\title{estimate parameters for COM-Poisson regression}

\description{ This package offers the ability to compute the COM-Poisson parameter estimates and associated standard errors.  This package also provides a hypothesis test for determining statistically significant data dispersion, and other model diagnostics.}

\details{
This package offers the ability to compute the COM-Poisson parameter estimates (via the cmp function) and associated standard errors.

Further, the user can perform a hypothesis test to determine the statistically significant need for using COM-Poisson regression to model the data.  The test addresses the matter of statistically significant dispersion.

The main order of functions is as follows:

1.  Compute Poisson estimates (using standard glm)
2.  Use Poisson estimates as starting values to determine COM-Poisson estimates (using cmp)
3.  Compute associated standard errors (using sdev function)

From here, there are lots of ways to proceed, so order doesn't matter:
 - Perform a hypothesis test to assess for statistically significant dispersion (using chisq and pval, or non-parametrically using parametric_bootstrap)
 - Compute leverage (using leverage) and deviance (using deviance)
 - Predict the outcome for new examples, using predict

}

\author{
Kimberly Sellers, Thomas Lotze
Maintainer: Thomas Lotze <thomas.lotze@thomaslotze.com>
}

\references{A Flexible Regression Model for Count Data, by Sellers & Shmueli, http://ssrn.com/abstract=1127359}

\examples{
# load freight data
  data(freight)

# Compute Standard Poisson estimates
  glm_model <- glm(broken ~ transfers, data=freight,
    family=poisson, na.action=na.exclude) # beta estimates
  print("The standard Poisson estimates for the beta vector are")
  print(coef(glm_model))

# Compute COM-Poisson estimates (under constant dispersion model)
start.time <- Sys.time()
  cmp_model = cmp(formula = broken ~ transfers, data=freight)
  print("The COM-Poisson estimates for the beta vector are")
  print(coef(cmp_model))
  print("The COM-Poisson estimate for the dispersion parameter nu is")
  print(nu(cmp_model))

# Compute associated standard errors for constant COM-Poisson estimates
  print("The associated standard errors for the betas in the constant dispersion case are")
  print(sdev(cmp_model))

# Perform likelihood ratio test for dispersion parameter
  # Test for dispersion equal or not equal to 1 (ie performing Poisson vs COM-Poisson regression)
  freight.chisq <- chisq(cmp_model)
  freight.pval <- pval(cmp_model)
  print(sprintf("The likelihood ratio chi-squared test statistic is \%0.5f
    and associated p-value (testing Poisson vs CMP regression) is \%0.5f",
    freight.chisq, freight.pval))

# Compute constant COM-Poisson leverage
  freight.lev <- leverage(cmp_model)
  print("The leverage of the points is")
  print(freight.lev)

# Compute constant COM-Poisson deviances
  # commented-out for speed
  # freight.CMPDev <- deviance(cmp_model)
  # print("The approximate constant dispersion standardized CMP Deviance is")
  # print(freight.CMPDev)

# Compute fitted values
  freight.fitted = predict(cmp_model, newdata=freight)
  print("The CMP fitted values are")
  print(freight.fitted)

# Compute residual values
  freight.constantCMPresids <- residuals(cmp_model)
  print("The CMP residuals are")
  print(freight.constantCMPresids)

# Compute MSE
  freight.constantCMP.MSE <- mean(freight.constantCMPresids^2)
  print("The MSE for the constant CMP regression is")
  print(freight.constantCMP.MSE)

# Compute predictions on new data
  new_data = data.frame(transfers=(0:10))
  freight.predicted = predict(cmp_model, newdata=new_data)
  plot(0:10, freight.predicted, type="l",
    xlab="number of transfers", ylab="predicted number broken")

# Compute parametric bootstrap results and use them to generate
#  0.95 confidence intervals for parameters
  #  commented-out for speed
  #  freight.CMPParamBoot <- parametric_bootstrap(cmp_model, n=1000)
  #  print(apply(freight.CMPParamBoot,2,quantile,c(0.025,0.975)))
}
