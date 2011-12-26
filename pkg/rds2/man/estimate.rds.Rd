\name{estimate.rds}
\alias{estimate.rds}
\title{Main function in rds2 package. Returns the ML estimate of population size and degree distribution.}
\usage{
  estimate.rds(data, init, const, arc = FALSE,
    maxit = 1000, ...)
}
\arguments{
  \item{func}{}

  \item{data}{}

  \item{init}{}

  \item{const}{}

  \item{arc}{}

  \item{maxit}{}

  \item{theta.range}{}

  \item{theta.minimum}{}

  \item{...}{}
}
\value{
  TBC
}
\description{
  Main function in rds2 package. Returns the ML estimate of
  population size and degree distribution.
}
\examples{
data(simulation)
estimate.rds(unlist(data3[1,]), const=20, arc=FALSE)
}
\author{
  johnros
}

