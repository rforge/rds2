\name{makeNjGrid}
\alias{makeNjGrid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Create grid of \eqn{\N[j]} starting values.
}
\description{
	Utility function for \code{\link{estimate.rds.grid}}. Creates a grid of \eqn{\N[j]} given a sample, \eqn{\theta} and \eqn{\beta} values.  
}
\usage{
makeNjGrid(sampled.degree.vector, theta.beta.grid, Nj.inflations)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sampled.degree.vector}{A numeric vector of the sampled degrees.}
  
  \item{theta.beta.grid}{The output of \code{\link{mekaBetaGrid}}}

  \item{Nj.inflations}{A numeric vector of inflation factors for each starting \eqn{\N[j]} vector.}
}
\value{
	Returns a list of starting values for \code{\link{estimates.rds.grid}}.
	}

\author{Jonathan D. Rosenblatt <john.ros.work@gmail.com>}

\examples{
Njs<- c(100,100,100,100); names(Njs)<- c("1","50","100","1000"); Njs<- as.table(Njs)
theta<- 1
beta<- 2.5e-8
tail(degree.sampled.vec<- generate.sample(theta, Njs, beta, sample.length=1e4))
theta.beta.grid<- makeBetaGrid(degree.sampled.vec, theta+c(-0.5, 0, 0.5), deflations= c(0.1, 0.5, 1))
makeNjGrid(degree.sampled.vec, theta.beta.grid, c(1,1000))
}

%\keyword{ ~kwd1 }

