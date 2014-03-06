\name{pval}
\alias{pval}

\title{P-value}
\description{a generic function for the p-value from the results of various model fitting functions. The function invokes particular methods which depend on the class of the first argument.} 

\usage{ pval(object, ...) }

\arguments{
  \item{object}{a model object}
  \item{...}{other parameters which might be required by the model}
}
\details{
See the documentation of the particular methods for details.
} 

\seealso{
       \code{\link{pval.cmp}}
}
\value{
The form of the value returned depends on the class of its argument. See the documentation of the particular methods for details of what is produced by that method.} 

\author{Thomas Lotze}