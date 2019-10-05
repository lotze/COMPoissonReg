#' Equidispersion Test
#' 
#' Likelihood ratio test for Equidispersion
#' 
#' @param object a model object
#' @param ... other parameters which might be required by the model
#' 
#' @details
#' A generic function for the likelihood ratio test for
#' equidispersion using the output of a fitted mode. The function invokes
#' particular methods which depend on the class of the first argument.
#' 
#' @return
#' Returns the test statistic and p-value determined from the \eqn{\chi.1^2}
#' distribution.
#' 
#' @author Thomas Lotze
#' @name equitest
#' @export
equitest = function(object, ...)
{
	UseMethod("equitest")
}

#' Leverage
#' 
#' A generic function for the leverage of points used in various model
#' fitting functions. The function invokes particular methods which
#' depend on the class of the first argument.
#' 
#' @param object a model object
#' @param ... other parameters which might be required by the model
#' 
#' @return
#' The form of the value returned depends on the class of its argument.
#' See the documentation of the particular methods for details of what is
#' produced by that method.
#' 
#' @details
#' See the documentation of the particular methods for details.
#' 
#' @author Thomas Lotze
#' @name leverage
#' @export
leverage = function(object, ...)
{
	UseMethod("leverage")
}

#' Estimate for dispersion parameter
#' 
#' A generic function for the dispersion parameter estimate from the results
#' of various model fitting functions. The function invokes particular methods
#' which depend on the class of the first argument.
#' 
#' @param object a model object
#' @param ... other parameters which might be required by the model
#' 
#' @details
#' See the documentation of the particular methods for details.
#' 
#' @return
#' The form of the value returned depends on the class of its argument. See
#' the documentation of the particular methods for details of what is
#' produced by that method.
#' 
#' @name nu
#' @export
nu = function(object, ...)
{
	UseMethod("nu")
}

#' Standard deviation
#' 
#' A generic function for standard deviation estimates from the results
#' of various model fitting functions. The function invokes particular
#' methods which depend on the class of the first argument.
#' 
#' @param object a model object
#' @param ... other parameters which might be required by the model
#' 
#' @return
#' The form of the value returned depends on the class of its argument.
#' See the documentation of the particular methods for details of what
#' is produced by that method.
#' 
#' @details
#' See the documentation of the particular methods for details.
#' 
#' @author Thomas Lotze
#' @name sdev
#' @export
sdev = function (object, ...)
{
	UseMethod("sdev")
}

#' Parametric Bootstrap
#' 
#' A generic function for the parametric bootstrap from the results of
#' various model fitting functions. The function invokes particular methods
#' which depend on the class of the first argument.
#' 
#' @param object a model object
#' @param ... other parameters which might be required by the model
#' @param reps Number of bootstrap repetitions.
#' @param report.period Report progress every \code{report.period} iterations.
#' 
#' @details
#' See the documentation of the particular methods for details.
#' 
#' @return
#' The form of the value returned depends on the class of its argument. See
#' the documentation of the particular methods for details of what is produced by that method.
#' 
#' @author Thomas Lotze
#' @name parametric.bootstrap
#' @export
parametric.bootstrap = function(object, reps = 1000, report.period = reps + 1, ...)
{
	UseMethod("parametric.bootstrap")
}
