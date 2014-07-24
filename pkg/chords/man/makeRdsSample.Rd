\name{makeRdsSample}
\alias{makeRdsSample}

\title{
Generate a synthetic (simulated) RDS sample.
}

\description{
Generates a sample from the sampling process assumed in the reference.
Well, actually, only the sufficient statistics required by \code{\link{estimate.b.k}} are returned.
}

\usage{
makeRdsSample(N.k, b.k, sample.length)
}

\arguments{
  \item{N.k}{An integer vector with the population frequency of each degree. }
  \item{b.k}{A numeric vector of the sampling rates of each degree. }
  \item{sample.length}{The length of the sample. Specified as the number of recruitees before termination.}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
  An object of class \code{rds-object} suitable for applying \code{\link{estimate.b.k}}.
}

\references{
  [1] Berchenko, Yakir, Jonathan Rosenblatt, and Simon D. W. Frost. "Modeling and Analyzing Respondent Driven Sampling as a Counting Process."" arXiv:1304.3505, April 11, 2013. HTTP://arXiv.org/abs/1304.3505.}

\author{
  Jonathan D. Rosenblatt \email{john.ros.work@gmail.com}
}

\note{
 The simulator does not prodice a whole RDS sample, but rather the sufficient statistics required for applying \code{\link{estimate.b.k}}.
}


\seealso{
\code{\link{estimate.b.k}}
}
\examples{
# setting up population data:
col.classes <- c('integer', 'integer', 'integer', 'integer', 'integer', 'integer', 'character')
f.loc <- paste(path.package('chords'), '/extdata/rds_data_Uganda02.csv',sep='')
rds.sample<- read.csv(file = f.loc, colClasses=col.classes)
rds.object<- initializeRdsObject(rds.sample)
rds.object$estimates <- estimate.b.k(rds.object = rds.object, impute.Nks = FALSE )

## Simulating rds.object:
scaler1 <- 10
rds.object$estimates$Nk.estimates <- rds.object$estimates$Nk.estimates* scaler1
rds.simulated.object <- makeRdsSample(
  N.k =rds.object$estimates$Nk.estimates , 
  b.k = exp(rds.object$estimates$log.bk.estimates),
  sample.length = 5e2)

# Estimating with simulated data:
rds.simulated.object$estimates <- estimate.b.k(rds.object = rds.simulated.object )
}
