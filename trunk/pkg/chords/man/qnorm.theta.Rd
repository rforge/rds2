\name{qnorm.theta}
\alias{qnorm.theta}
\alias{inv.qnorm.theta}
\title{ Maps \eqn{\theta} to the real line and vise versa.}
\usage{
  qnorm.theta(theta, ... )
    
   inv.qnorm.theta(qnorm.theta, ...)
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
  Jonathan D. Rosenblatt \email{john.ros@gmail.com}
}
\seealso{
  \code{\link{estimate.rds}}  
}

