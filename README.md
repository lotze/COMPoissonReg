COMPoissonReg
=============
The COMPoissonReg R package fits COM-Poisson regression (Sellers & Shmueli, 2010) and zero-inflated (ZI) COM-Poisson regression models (Sellers & Raim, 2016). The syntax to fit a model looks like
```R
cmp.out <- glm.cmp(formula.lambda = y ~ x1 + x2, formula.nu = ~ x1, formula.p = ~ x2)
```
where `formula.lambda` is a regression model for the COM-Poisson rate parameter, `formula.nu` is a model for the dispersion parameter, and `formula.p` is a model for the zero-inflation parameter. See the package documentation for details.

## References
* Kimberly F. Sellers & Galit Shmueli (2010). A Flexible Regression Model for
Count Data. Annals of Applied Statistics, 4(2), 943-961. [[http]](http://projecteuclid.org/euclid.aoas/1280842147)
* Kimberly F. Sellers and Andrew M. Raim (2016). A Flexible Zero-Inflated Model
to Address Data Dispersion. Computational Statistics and Data Analysis, 99,
68-80. [[http]](http://www.sciencedirect.com/science/article/pii/S0167947316000165)