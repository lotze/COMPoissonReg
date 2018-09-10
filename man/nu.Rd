\name{nu}
\alias{nu}

\title{Estimate for dispersion parameter}
\description{a generic function for the dispersion parameter estimate from the results of various model fitting functions. The function invokes particular methods which depend on the class of the first argument.} 

\usage{ nu(object, ...) }

\arguments{
  \item{object}{a model object}
  \item{...}{other parameters which might be required by the model}
}
\details{
See the documentation of the particular methods for details.
} 

\seealso{
       \code{\link{nu.cmp}},
       \code{\link{nu.zicmp}}
}
\value{
The form of the value returned depends on the class of its argument. See the documentation of the particular methods for details of what is produced by that method.} 

\author{Thomas Lotze}