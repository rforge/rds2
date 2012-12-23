\name{transform.theta}
\alias{transform.theta}
\alias{inv.transform.theta}
\title{ Maps \eqn{\theta} to the real line and vise versa.}
\usage{
  transform.theta(theta, ... )
    
   inv.transform.theta(qnorm.theta, ...)
}
\arguments{
  \item{theta}{\eqn{\theta } in its original scale}
  
  \item{qnorm.theta}{\eqn{\theta } rescaled by \code{qnorm.theta}}
  
  \item{...}{ For future usage.}

  
}
\value{
  A rescaled \eqn{\theta}.
}
\description{
  Maps \eqn{\theta} to the real line and vise versa.
}
\note{
These functions were created in the view that theta values might be constrained is some way. Eventually, this was not done so at present they are both identity functions.
}

\author{
  Jonathan D. Rosenblatt \email{john.ros.work@gmail.com}
}
\seealso{
  \code{\link{estimate.rds}}  
}

