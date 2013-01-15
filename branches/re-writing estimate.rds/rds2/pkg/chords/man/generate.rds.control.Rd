\name{generate.rds.control}

\alias{generate.rds.control}

\title{Generator for parameters controlling \code{\link{estimate.rds}}}
\description{
Generator for parameters controlling \code{\link{estimate.rds}}
}
\usage{
generate.rds.control(maxit = 500, method = "Nelder-Mead", Nj.inflations = c(1, 2, 10), beta.inflations = exp(-seq(2, 0, by = -0.02)), fnscale = -1000, beta.scale = 1e+06)
}

\arguments{
  \item{maxit}{
Number of maximal optimization iterations passed to \code{\link{optim}}.
}
  \item{method}{
Optimization methods passed to \code{\link{optim}}.
}
  \item{Nj.inflations}{
Inflation factors for creating set of \eqn{N[j]} initialization values. 
}
  \item{beta.inflations}{
Inflation factors for creating set of \eqn{\beta} initialization values.
}
  \item{fnscale}{
Target function scale. Passed to \code{\link{optim}}.
}
  \item{beta.scale}{
\eqn{\beta} parameter scale. 
}
}
\value{
Wraps its arguments as a named list needed for \code{\link{estimate.rds}}}.
\author{
Jonathan Rosenblatt<john.ros.work@gmail.com>
}

\examples{
generate.rds.control()
}
