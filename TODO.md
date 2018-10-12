# To Do

We currently don't throw errors if the z-function evaluates to infinity.
This causes some weird results, like the parameters lambda = exp(5.25),
nu = 0.4 always causing rcmp to draw 203 as the value.

Add ymax as an option (not global) to regression functions

