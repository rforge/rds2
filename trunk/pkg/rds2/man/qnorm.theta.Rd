\name{qnorm.theta}
\alias{qnorm.theta}
\alias{inv.qnorm.theta}
\title{ Maps \eqn{\theta} to the real line and vise versa.}
\usage{
  qnorm.theta(theta, const, range = 20, minimum = -10)
    
   inv.qnorm.theta(qnorm.theta, const, range = 20,
    minimum = -10)
}
\arguments{
  \item{theta}{\eqn{\theta } in its original scale}
  
  \item{qnorm.theta}{\eqn{\theta } rescaled by \code{qnorm.theta}}

  \item{const}{Control the concentration around zero of the mapping of the theta parameter from the constrained domain to the real line. Typically passed by \code{estimate.rds}}

 \item{range}{ Defines the range of \eqn{\theta} from \code{minimum} to \code{minimum+range} }
 
 \item{minimum}{ The minimal value for \eqn{\theta} in the original scale.} 
  
}
\value{
  A rescaled \eqn{\theta}.
}
\description{
  Maps \eqn{\theta} to the real line and vise versa.
}
\author{
  Jonathan D. Rosenblatt \email{john.ros@gmail.com}
}
\seealso{
  \code{\link{estimate.rds}}  
}

