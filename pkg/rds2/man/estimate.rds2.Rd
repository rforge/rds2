\name{estimate.rds2}
\alias{estimate.rds2}
\title{Same as estimate.rds, but for a fixed theta.}
\usage{
  estimate.rds2(data, Sij, init, const, arc = FALSE,
    maxit = 10000, theta, ...)
}
\arguments{
  \item{data}{}

  \item{Sij}{}

  \item{init}{}

  \item{const}{}

  \item{arc}{}

  \item{maxit}{}

  \item{theta}{}

  \item{...}{}
}
\value{
  TBC
}
\description{
  Same as estimate.rds, but for a fixed theta.
}
\examples{
data(brazil)
estimate.rds2(data=data.degree, , Sij = data.Sjt, const=50, arc=FALSE, maxit=1000, theta = 1)
}
\author{
  johnros
}

