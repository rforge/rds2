\name{inv.qnorm.theta}
\alias{inv.qnorm.theta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
inv.qnorm.theta(qnorm.theta, const, range = 20, minimum = -10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{qnorm.theta}{
%%     ~~Describe \code{qnorm.theta} here~~
}
  \item{const}{
%%     ~~Describe \code{const} here~~
}
  \item{range}{
%%     ~~Describe \code{range} here~~
}
  \item{minimum}{
%%     ~~Describe \code{minimum} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
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
function (qnorm.theta, const, range = 20, minimum = -10) 
{
    normalized.theta <- pnorm(qnorm.theta/const)
    normalized.theta2 <- (normalized.theta * range) + minimum
    return(normalized.theta2)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
