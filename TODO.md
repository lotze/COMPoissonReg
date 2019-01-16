# To Do

We currently don't throw errors if the z-function evaluates to infinity.
This causes some weird results, like the parameters lambda = exp(5.25),
nu = 0.4 always causing rcmp to draw 203 as the value.

Add ymax as an option (not global) to regression functions

There is an issue with formula processing. Here is a small example to
illustrate.
``` {r}
y <- rcmp(250, lambda = 10, nu = 0.95)
# Doesn't work
out <- glm.cmp(y ~ 1)
out <- glm.cmp(y ~ 1, formula.nu = ~ 1)
# Workaround
out <- glm.cmp(y ~ 1, formula.nu = y ~ 1)
```

