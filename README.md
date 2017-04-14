COMPoissonReg
=============
The COMPoissonReg R package fits COM-Poisson regression (Sellers & Shmueli, 2010) and zero-inflated (ZI) COM-Poisson regression models (Sellers & Raim, 2016). The general ZI COM-Poisson regression model is
$$
Y_i \sim p_i I(y_i = 0) + (1 - p_i) \frac{\lambda_i^y}{(y!)^{\nu_i} Z(\lambda_i, \nu_i)}, \quad y_i = 0, 1, \ldots, \\
\log \lambda_i = x_i^T \beta, \quad
\log \nu_i = s_i^T \gamma, \quad
\log\left[ \frac{p_i}{1-p_i} \right] = w_i^T \zeta,
$$
where $Y_i$ are count variables, $(x_i, s_i, w_i)$ are predictors, and $(\beta, \gamma, \zeta)$ are regression coefficients. The syntax to fit a model looks like this.
```R
cmp.out <- glm.cmp(y ~ x1 + x2, formula.nu = ~ x1, formula.p = ~ x2)
```
See the package documentation for details.

## References
* Kimberly F. Sellers & Galit Shmueli (2010). A Flexible Regression Model for
Count Data. Annals of Applied Statistics, 4(2), 943-961. [[http]](http://projecteuclid.org/euclid.aoas/1280842147)
* Kimberly F. Sellers and Andrew M. Raim (2016). A Flexible Zero-Inflated Model
to Address Data Dispersion. Computational Statistics and Data Analysis, 99,
68-80. [[http]](http://www.sciencedirect.com/science/article/pii/S0167947316000165)