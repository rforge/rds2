\name{estimate.rds}
\alias{estimate.rds}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
estimate.rds(data, init, const, arc, maxit = 20000, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{init}{
%%     ~~Describe \code{init} here~~
}
  \item{const}{
%%     ~~Describe \code{const} here~~
}
  \item{arc}{
%%     ~~Describe \code{arc} here~~
}
  \item{maxit}{
%%     ~~Describe \code{maxit} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (data, init, const, arc, maxit = 20000, ...) 
{
    N.j <- rep(0, max(data))
    uniques <- unique(data)
    uniques <- uniques[uniques != 0]
    s.uniques <- sort(uniques)
    param.size <- length(uniques)
    data.table <- table(data)[-1]
    likelihood.wrap1 <- function(par) {
        final.result <- -Inf
        ccc <- exp(par[1])
        theta <- inv.qnorm.theta(par[2], const = const)
        N.j <- rep(0, max(data))
        N.j[s.uniques] <- exp(tail(par, -2))
        if (any(N.j[s.uniques] < data.table)) 
            return(final.result)
        if (is.numeric(ccc) && is.numeric(theta) && !is.infinite(theta) && 
            !is.infinite(ccc)) {
            result <- .C("likelihood", sample = as.integer(data), 
                c = as.numeric(ccc), theta = as.numeric(theta), 
                Nj = as.numeric(N.j), constant = as.numeric(const), 
                n = as.integer(length(data)), N = as.integer(length(N.j)), 
                arc = arc, result = as.double(0))
            final.result <- result$result
        }
        return(final.result)
    }
    cap <- max(data)
    cap.arcs <- sum.ranks(data)
    max.cap.arcs <- max(cap.arcs)
    log.c <- -log(max(sum.ranks(data) * sum(data)) * max(data)) - 
        1
    if (missing(init)) {
        init <- list(six0.2 = c(log.c = log.c, qnorm.theta = qnorm.theta(0.2, 
            const), log(data.table) + 1), six1 = c(log.c = log.c, 
            qnorm.theta = qnorm.theta(1, const), log(data.table) + 
                1), six2 = c(log.c = log.c, qnorm.theta = qnorm.theta(2, 
            const), log(data.table) + 1))
    }
    likelihood.optim <- lapply(init, function(x) {
        try(optim(par = x, fn = likelihood.wrap1, control = list(fnscale = -1, 
            maxit = maxit), ...))
    })
    prepare.result <- function(x) {
        N.j[s.uniques] <- exp(tail(x$par, -2))
        result <- list(c = exp(x$par[[1]]), theta = inv.qnorm.theta(x$par[[2]], 
            const), Nj = N.j, x)
        return(result)
    }
    return(lapply(likelihood.optim, function(x) try(prepare.result(x))))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
