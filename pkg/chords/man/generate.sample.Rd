\name{generate.sample}
\alias{generate.sample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generates a random respondent driven sample  
}
\description{
Generates a discrete-time respondent driven sample assuming the model in the reference. Typically used to test \code{\link{estimate.rds}}
}
\usage{
generate.sample(theta, Njs, beta, sample.length, double.sampling.warning.thresh = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
See reference.
}
  \item{Njs}{
See reference.
}
  \item{beta}{
See reference.
}
  \item{sample.length}{
The duration of sampling in discrete time units.
}
  \item{double.sampling.warning.thresh}{
Warns when the parameters are such that the probability of double sampling in a single time step is non negligible.}
}
\value{
A numeric vector of the degree of the individual sampled at each time step.}
\references{
	Respondent Driven Sampling as an Epidemic Process
	Berchenko Y., Rosenblatt J.D., White R.G.,Frost S.D.W. (2012) 
  \emph{?}

}
\author{
  Jonathan D. Rosenblatt \email{john.ros.work@gmail.com}
}


\seealso{
\code{\link{estimate.rds}}
}
\examples{
network.size<- 2000
network.density<- 0.002
neighbour.count<- lapply(seq(1, network.size), function(x) rbinom(1, network.size, network.density))

Njs<- table(unlist(neighbour.count[neighbour.count>0]))
theta<- 1
beta<- 9e-8 
degree.sampled.vec<- generate.sample(theta = theta, Njs = Njs, beta = beta, sample.length = 5000)
plot(degree.sampled.vec, type='h')
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datagen}
\keyword{survey}
