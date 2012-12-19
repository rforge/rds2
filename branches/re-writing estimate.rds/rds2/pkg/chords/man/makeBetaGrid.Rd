\name{makeBetaGrid}
\alias{makeBetaGrid}
\title{
 Create grid of \eqn{\beta} starting values.
}
\description{
	Utility function for \code{\link{estimate.rds.grid}}. Creates a grid of \eqn{\beta} given a sample, \eqn{\theta} and \eqn{\beta} values.  
}
\usage{
makeBetaGrid(sampled.degree.vector, thetas, betas, beta.inflations)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sampled.degree.vector}{A numeric vector of the sampled degrees.}
  
  \item{thetas}{A numeric vector of starting \eqn{\theta} values.}
  
  \item{betas}{A numeric vector of \eqn{\beta} values to be deflated and returned in the needed format for \code{\link{estimate.rds.grid}}}
  
  \item{beta.inflations}{A numeric vector of deflation factors for each starting \eqn{\beta}.}
}
\value{
	Returns a list of starting values for \code{\link{estimate.rds.grid}}.
}
\author{
Jonathan D. Rosenblatt <john.ros.work@gmail.com>
}

\examples{
Njs<- c(100,100,100,100); names(Njs)<- c("1","50","100","1000"); Njs<- as.table(Njs)
theta<- 1
beta<- 2.5e-8
tail(degree.sampled.vec<- generate.sample(theta, Njs, beta, sample.length=1e4))
makeBetaGrid(degree.sampled.vec, theta+c(-0.5, 0, 0.5), deflations= c(0.1, 0.5, 1))
}

%\keyword{ ~kwd1 }

