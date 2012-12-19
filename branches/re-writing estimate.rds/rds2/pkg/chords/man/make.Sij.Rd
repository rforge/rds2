\name{make.Sij}
\alias{make.Sij}
\title{Make the \eqn{Sij} matrix assuming no drop outs during sampling }
\usage{
  make.Sij(data)
}
\arguments{
  \item{data}{ An integer vector of the degrees of each individual sampled.}
}
\value{
  A matrix with the number of individuals with degree i (as row names) at sampling stage j (as column). 
}
\description{
  Make the \eqn{Sij} matrix from the sample path, assuming no drop outs during sampling.
}
\author{
  Jonathan D. Rosenblatt \email{john.ros.work@gmail.com}
}
\seealso{
  \code{\link{estimate.rds}}  
}

