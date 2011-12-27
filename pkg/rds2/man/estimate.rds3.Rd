\name{estimate.rds3}
\alias{estimate.rds3}
\title{Main function in rds2 package. Returns the ML estimate of population size and degree distribution.}
\usage{
  estimate.rds3(data, Sij, init, const, arc = FALSE,
    maxit = 10000, initial.thetas, theta.minimum,
    theta.range)
}
\arguments{
  \item{data}{ An integer vector of the degrees of each individual sampled.}

  \item{Sij}{ An integer matrix of the counts of individuals with rank j at sampling period t. Row names correspond to the degrees.}

  \item{init}{Optional initialization values.}

  \item{const}{Control the concentration of the mapping of the theta parameter from the contstrained domain to the real line.}

  \item{arc}{Should the degrees in the snowball be taken in account (TRUE) or just the size of the snowball (FALSE). Defaults to FALSE.}

  \item{maxit}{Number of maximal optimization iterations.}

  \item{initial.thetas}{ Vector of initial values of theta parameter.}

  \item{theta.minimum}{Lower bound on the allowed values of theta.}

  \item{theta.range}{Allowed range of theta values. From theta.minimum to theta.minimum+range.}

}
\description{
  Performs maximum likelihood estimation of the population
  size, degree distribution and theta [explain] allowing
  (a) "soft" constraints on theta and (b) an input Sij
  matrix.
}
\examples{
data(simulation)
temp.data<- unlist(data3[1,7000:7500])
estimate.rds3(temp.data, Sij = make.Sij(temp.data), initial.thetas = c(1,10), arc = FALSE, maxit = 1000, const = 0.5, theta.minimum = -0.5, theta.range = 2)
}
\author{
  Jonathan Rosenblatt
}

