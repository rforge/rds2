\name{estimate.rds2}
\alias{estimate.rds2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
estimate.rds2(data, Sij, init, const, arc, maxit = 10000, theta, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{Sij}{
%%     ~~Describe \code{Sij} here~~
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
  \item{theta}{
%%     ~~Describe \code{theta} here~~
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
function (data, Sij, init, const, arc, maxit = 10000, theta, 
    ...) 
{
    N.j <- rep(0, max(data))
    uniques <- unique(data)
    uniques <- uniques[uniques != 0]
    s.uniques <- sort(uniques)
    param.size <- length(uniques)
    data.table <- table(data)[-1]
    S <- compute.S(Sij)
    likelihood.wrap1 <- function(par) {
        final.result <- -Inf
        ccc <- exp(par[1])
        N.j <- rep(0, max(data))
        N.j[s.uniques] <- exp(tail(par, -1))
        if (any(N.j[s.uniques] < data.table)) 
            return(final.result)
        if (is.numeric(ccc) && is.numeric(theta) && !is.infinite(theta) && 
            !is.infinite(ccc)) {
            result <- .C("likelihood", sample = as.integer(data), 
                Sij = as.integer(as.matrix(Sij)), S = as.integer(S), 
                c = as.numeric(ccc), theta = as.numeric(theta), 
                Nj = as.numeric(N.j), constant = as.numeric(const), 
                observed_degrees = as.integer(rownames(Sij)), 
                n = as.integer(length(data)), N = as.integer(length(N.j)), 
                N_observed = as.integer(nrow(Sij)), arc = arc, 
                result = as.double(0))
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
        init <- list(six = c(log.c = log.c, log(data.table) + 
            1))
    }
    likelihood.optim <- lapply(init, function(x) {
        try(optim(par = x, fn = likelihood.wrap1, control = list(fnscale = -1, 
            maxit = maxit), ...))
    })
    prepare.result <- function(x) {
        N.j[s.uniques] <- exp(tail(x$par, -1))
        result <- list(c = exp(x$par[[1]]), theta = theta, Nj = N.j, 
            x)
        return(result)
    }
    return(lapply(likelihood.optim, function(x) try(prepare.result(x))))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
