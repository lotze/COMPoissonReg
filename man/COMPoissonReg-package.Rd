\name{COMPoissonReg-package}

\alias{COM-PoissonReg-package}
\alias{COM-PoissonReg}
\alias{COMPoissonReg}

\docType{package}

\title{Estimate parameters for COM-Poisson regression}

\description{ This package offers the ability to compute the parameter estimates for a COM-Poisson or zero-inflated (ZI) COM-Poisson regression and associated standard errors.  This package also provides a hypothesis test for determining statistically significant data dispersion, and other model diagnostics.}

\details{
This package offers the ability to compute the COM-Poisson parameter estimates and associated standard errors for a regular regression model or a zero-inflated regression model (via the glm.cmp function).

Further, the user can perform a hypothesis test to determine the statistically significant need for using COM-Poisson regression to model the data.  The test addresses the matter of statistically significant dispersion.

The main order of functions for COM-Poisson regression is as follows:
\enumerate{
\item  Compute Poisson estimates (using glm for Poisson regression or pscl for
ZIP regression)
\item  Use Poisson estimates as starting values to determine COM-Poisson estimates (using glm.cmp)
\item  Compute associated standard errors (using sdev function)
}

From here, there are lots of ways to proceed, so order is irrelevant:
\itemize{
\item Perform a hypothesis test to assess for statistically significant dispersion (using equitest or parametric_bootstrap)
\item Compute leverage (using leverage) and deviance (using deviance)
\item Predict the outcome for new examples, using predict
}

The package also supports fitting of the zero-inflated COM-Poisson model (ZICMP).
Most of the tools available for COM-Poisson are also available for ZICMP.
}

\author{
Kimberly Sellers, Thomas Lotze (Maintainer, <thomas.lotze@thomaslotze.com>),
Andrew M. Raim
}

\references{
Kimberly F. Sellers & Galit Shmueli (2010). A Flexible Regression Model for
Count Data. Annals of Applied Statistics, 4(2), 943-961.

Kimberly F. Sellers and Andrew M. Raim (2016). A Flexible Zero-Inflated Model
to Address Data Dispersion. Computational Statistics and Data Analysis, 99,
68-80.
}

\examples{
  ## load freight data
  data(freight)

  # Fit standard Poisson model
  glm.out <- glm(broken ~ transfers, data=freight,
    family=poisson, na.action=na.exclude)
  print(glm.out)

  # Fit COM-Poisson model (with intercept-only regression linked to the
  # dispersion parameter)
  cmp.out <- glm.cmp(broken ~ transfers, data=freight)
  print(cmp.out)
  coef(cmp.out)
  nu(cmp.out)[1]

  # Compute associated standard errors
  sdev(cmp.out)

  # Get the full covariance matrix for the estimates
  vcov(cmp.out)

  # Likelihood ratio test for dispersion parameter
  # Test for H_0: dispersion equal to 1 vs. H_1: not equal to 1
  # (i.e. Poisson vs. COM-Poisson regression models)
  lrt <- equitest(cmp.out)

  # Compute constant COM-Poisson leverage
  lev <- leverage(cmp.out)

  \donttest{
  # Compute constant COM-Poisson deviances
  dev <- deviance(cmp.out)
  }

  # Compute fitted values
  y.hat <- predict(cmp.out, newdata=freight)

  # Compute residual values
  res <- residuals(cmp.out)
  print(summary(res))

  # Compute MSE
  mean(res^2)

  # Compute predictions on new data
  new_data <- data.frame(transfers=(0:10))
  y.hat <- predict(cmp.out, newdata=new_data)
  plot(0:10, y.hat, type="l",
    xlab="number of transfers", ylab="predicted number broken")

  \donttest{
  # Compute parametric bootstrap results and use them to generate
  # 0.95 confidence intervals for parameters.
  cmp.boot <- parametric_bootstrap(cmp.out, reps=1000)
  print(apply(cmp.boot, 2, quantile, c(0.025,0.975)))
  }

  \donttest{
  ## load couple data
  data(couple)

  # Fit standard Poisson model
  glm.out <- glm(UPB ~ EDUCATION + ANXIETY, data=couple, family=poisson)
  print(glm.out)

  # Fit ZICMP model
  zicmp.out <- glm.cmp(UPB ~ EDUCATION + ANXIETY,
    formula.nu = ~ 1,
    formula.p = ~ EDUCATION + ANXIETY,
    data=couple)
  print(zicmp.out)

  # Compute standard errors for estimates of coefficients
  sdev(zicmp.out)

  # Get the full covariance matrix for the estimates
  vcov(zicmp.out)

  # Likelihood ratio test for equidispersion (H0: nu = 1 vs H1: not)
  equitest(zicmp.out)

  # Compute fitted values
  y.hat <- predict(zicmp.out)

  # Compute residuals
  res.raw <- residuals(zicmp.out, type = "raw")
  res.quan <- residuals(zicmp.out, type = "quantile")
  print(summary(res.raw))
  print(summary(res.quan))

  # Compute predictions on new data
  new_data <- data.frame(EDUCATION = round(1:20 / 20), ANXIETY = seq(-3,3, length.out = 20))
  y.hat.new <- predict(zicmp.out, newdata=new_data)
  print(y.hat.new)
  
  # Compute parametric bootstrap results and use them to generate
  # 0.95 confidence intervals for parameters.
  zicmp.boot <- parametric_bootstrap(zicmp.out, reps=1000)
  print(apply(zicmp.boot, 2, quantile, c(0.025,0.975)))
  }
}
