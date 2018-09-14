# To Do

Improving computation of z-function
* Both z_approx and z_exact functions give infinite logz values when nu is
  very close to zero. We weren't seeing this before, since we truncated the
  sum in the z function.
* Should we add back functions for truncated versions of CMP?
* We currently don't throw errors if the z-function evaluates to infinity.
  This causes some weird results, like the parameters lambda = exp(5.25),
  nu = 0.4 always causing rcmp to draw 203 as the value.
* For couple example, FIM and Hessian seem to give very different vcovs ...
* Retest everything
* Why is equitest's teststat = 0 for the couple example? Is that right?

Add documentation for global options

Add documentation for mc.reps in vcov and sdev, if we'll keep it

Problem with rcmp
lambda <- 10
nu <- 0.2
