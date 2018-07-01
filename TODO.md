# To Do

Improving computation of z-function
* May need adaptive code to compute other z-like functions.
* Or... put the max argument back in
* Other z-like functions may need very large max for accuracy. On the other
  hand, taking the numerical derivative is complicated by the fact that their
  parameters are bounded at zero. One option is to compute the derivatives
  using a transformation to Euclidean space (but this still might lead to
  problems).
* For couple example, FIM and Hessian seem to give very different vcovs ...
* May not want to use FIM as default variance method anymore (consider Hessian)
* Retest everything
* May want to delete C++ cmp_expected_value - numDeriv seems to be more stable

