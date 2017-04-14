\name{parametric_bootstrap}
\alias{parametric_bootstrap}

\title{Parametric Bootstrap}
\description{a generic function for the parametric bootstrap from the results of various model fitting functions. The function invokes particular methods which depend on the class of the first argument.} 

\usage{ parametric_bootstrap(object, reps = 1000, report.period = reps+1, ...) }

\arguments{
  \item{object}{a model object}
  \item{...}{other parameters which might be required by the model}
  \item{reps}{Number of bootstrap repetitions.}
  \item{report.period}{Report progress every \code{report.period} iterations.}
}
\details{
See the documentation of the particular methods for details.
} 

\seealso{
       \code{\link{parametric_bootstrap.cmp}},
       \code{\link{parametric_bootstrap.zicmp}}
}
\value{
The form of the value returned depends on the class of its argument. See the documentation of the particular methods for details of what is produced by that method.} 

\author{Thomas Lotze}
