\name{compute.S}
\alias{compute.S}
\title{ A snow ball constructor}
\usage{
  compute.S(Sij)
}
\arguments{
  \item{Sij}{ An integer matrix of the counts of individuals with rank i at sampling period j. Row names correspond to the degrees.}
}
\value{
  \item{S}{A numeric vector of the cumulative number of individuals sampled}
}
\description{
  Constructs the cumulative number of individuals sampled (\code{S}) given the 
}
\author{
  Jonathan D. Rosenblatt \email{john.ros@gmail.com}
}
\seealso{
  \code{\link{estimate.rds}}  
}

