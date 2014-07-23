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

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (N.k, b.k, sample.length) 
{
    n.ks <- rep(0, length(N.k))
    event.time <- 0
    I.t <- 1
    rds.object.simulated <- rdsObjectConstructor(rds.sample = data.frame(NS1 = rep(NA, 
        sample.length), interviewDt = rep(NA, sample.length)), 
        I.t = rep(NA, sample.length), degree.in = rep(0, sample.length), 
        degree.out = rep(0, sample.length))
    for (period in seq_len(sample.length)) {
        lambda.k <- updateLambdas(N.k, n.ks, b.k, I.t)
        time.to.event <- rexp(1, sum(lambda.k))
        event.time <- event.time + time.to.event
        rds.object.simulated$rds.sample$interviewDt[period] <- event.time
        event.type <- rmultinom(1, 1, prob = lambda.k)
        in.degree <- which(as.logical(event.type))
        rds.object.simulated$degree.in[period] <- in.degree
        n.ks <- n.ks + event.type
        rds.object.simulated$I.t[period] <- I.t
        I.t <- I.t + 1
        rds.object.simulated$rds.sample$NS1[period] <- in.degree
    }
    return(rds.object.simulated)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
