\name{estimate.rds2}
\alias{estimate.rds2}
\title{Same as estimate.rds, but for a fixed theta.}
\usage{
  estimate.rds2(data, Sij, init, const, arc = FALSE,
    maxit = 10000, theta)
}
\arguments{
  \item{data}{ An integer vector of the degrees of each individual sampled.}

  \item{Sij}{ An integer matrix of the counts of individuals with rank j at sampling period t. Row names correspond to the degrees.}

  \item{init}{Optional initialization values.}

  \item{const}{Control the concentration of the mapping of the theta parameter from the contstrained domain to the real line.}

  \item{arc}{Should the degrees in the snowball be taken in account (TRUE) or just the size of the snowball (FALSE). Defaults to FALSE.}

  \item{maxit}{Number of maximal optimization iterations.}

  \item{theta}{Value of theta parameter. Assumed known.}
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

