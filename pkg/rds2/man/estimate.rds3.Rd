\name{estimate.rds3}
\alias{estimate.rds3}
\title{Main function in rds2 package. Returns the ML estimate of population size and degree distribution.}
\usage{
  estimate.rds3(data, Sij, init, const, arc = FALSE,
    maxit = 10000, initial.thetas, theta.minimum,
    theta.range, ...)
}
\arguments{
  \item{data}{}

  \item{Sij}{}

  \item{init}{}

  \item{const}{}

  \item{arc}{}

  \item{maxit}{}

  \item{initial.thetas}{}

  \item{theta.minimum}{}

  \item{theta.range}{}

  \item{...}{}
}
\description{
  Performs maximum likelihood estimation of the population
  size, degree distribution and theta [explain] allowing
  (a) "soft" constraints on theta and (b) an input Sij
  matrix.
}
\author{
  Jonathan Rosenblatt
}

