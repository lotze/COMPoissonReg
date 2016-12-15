\name{equitest}
\alias{equitest}

\title{Likelihood ratio test for Equidispersion}
\description{A generic function for the likelihood ratio test for
equidispersion using the output of a fitted mode. The function invokes
particular methods which depend on the class of the first argument.}

\usage{ equitest(object, ...) }

\arguments{
  \item{object}{a model object}
  \item{...}{other parameters which might be required by the model}
}
\details{
See the documentation of the particular methods for details.
}

\seealso{
       \code{\link{equitest.cmp}}
       \code{\link{equitest.zicmp}}
}
\value{
The form of the value returned depends on the class of its argument. See the documentation of the particular methods for details of what is produced by that method.}

\author{Thomas Lotze}
