\name{COMPoissonReg-options}
\alias{COMPoissonReg-options}

\title{Package options}
\description{Global options used by the COMPoissonReg package.} 

\details{
	\code{options(COMPoissonReg.optim.method = 'L-BFGS-B')}
	\code{options(COMPoissonReg.optim.control = list(maxit = 150))}
	\code{options(COMPoissonReg.grad.eps = 1e-5)}
	\code{options(COMPoissonReg.hess.eps = 1e-2)}
}

\arguments{
  \item{COMPoissonReg.optim.method}{Optim method to use when computing
  maximum likelihood estimates.}
  \item{COMPoissonReg.optim.control}{A list to be passed to \code{control}
  when calling \code{optim}. \code{fnscale} will be ignored if specified.}
  \item{COMPoissonReg.grad.eps}{Distance to be used when finite differences
  are taken.}
  \item{COMPoissonReg.hess.eps}{Distance to be used when finite second
  differences are taken.}
}
