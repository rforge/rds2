\name{initializeRdsObject}
\alias{initializeRdsObject}

\title{
Construct a \code{rds-object} from a data.frame.
}
\description{
Given a data frame with the appropriate variables, initializes a \code{rds-object} with the components required by the \code{\link{estimate.b.k}} function for estimation.
}
\usage{
initializeRdsObject(rds.sample)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rds.sample}{A data frame with required columns. See Details. }
}
\details{
The essence of the function is in recovering the sampling snowball required by \code{\link{estimate.b.k}}. 
The function allows for recruiters to enter and exit the sampling snowball.
It currently assumes there is a single seed recruiter (the snowball at time 0). This is because the "vanilla" RDS file format does not include recruiter information.

The \code{rds.sample} object is assumed to be a data frame.
Following the RDS file specification in [TODO: reference], the following column names are assumed:
\enumerate{
\item{MyUniID}{an identifier of the sampling unit.[not required]}
\item{NS1}{The reported degree.[required]}
\item{refCoupNum}{ The number of the referring coupon.}
\item{coup1}{The number of the 1st supplied coupon. NA if none. [required].}
\item{coupX}{The number of the Xth supplied coupon. NA if none.[not required]}
\item{interviewDt}{The time of the interview. In numeric representation from some origin. Ties are not allowed.}
}
}

\value{
A list with the following components.
\itemize{
\item{rds.sample}{The input data frame. After ordering along time of arrival.}
\item{I.t}{The sampling snowball. A list including the following items:
\item{I.t}{An integer of the count of the sampling individuals at the moments of recruitment.}
\item{degree.in}{An integer with the degree of an added recruiter at the moments of recruitment.}
\item{degree.out}{An integer with the degree of a removed recruiter at the moment of recruitment.}}
\item{original.ordering}{The order of the arrivals as was inputed in \code{rds.sample$interviewDt}}
\item{estimates}{A placeholder for the future output of \code{\link{estimate.b.k}}}
}
}

\references{
  [1] Berchenko, Yakir, Jonathan Rosenblatt, and Simon D. W. Frost. "Modeling and Analyzing Respondent Driven Sampling as a Counting Process."" arXiv:1304.3505, April 11, 2013. HTTP://arXiv.org/abs/1304.3505.
}
\author{
  Jonathan D. Rosenblatt \email{john.ros.work@gmail.com}
}

\seealso{
  \code{\link{estimate.b.k}}, \code{\link{makeRdsSample}}
}
\examples{
# See examples in estimate.b.k()
}