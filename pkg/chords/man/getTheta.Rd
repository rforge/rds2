\name{getTheta}
\alias{getTheta}

\title{
Recover the "discoverability coefficient".
}
\description{
Estimates the effect of the degree on the rate of sampling. Also known as the "coefficient of discoverability" in the oil-discovery literature [2].
Formally, we estimate \eqn{theta} and \eqn{beta_0} assuming that \eqn{beta_k:= beta_0 * k^theta}.
}

\usage{
  getTheta(rds.object, robust=TRUE)
}

\arguments{
  \item{rds.object}{A \code{rds-object} with a \code{estimates} component as returned by \code{\link{estimate.b.k}} }
  \item{robust}{Should \eqn{beta_0} and \eqn{theta} be recovered from \eqn{beta_k} using a robust method (default) or not.}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
A list including the following components:
\item{log.beta_0}{The log of \eqn{beta_0} in \eqn{beta_k:= beta_0 * k^theta}.}
\item{theta}{\eqn{theta} in \eqn{beta_k:= beta_0 * k^theta}.}
}

\references{
  [1] Berchenko, Yakir, Jonathan Rosenblatt, and Simon D. W. Frost. "Modeling and Analyzing Respondent Driven Sampling as a Counting Process."" arXiv:1304.3505. HTTP://arXiv.org/abs/1304.3505.
  [2] Bloomfield, P., K.S. Deffeyes, B. Silverman, G.S. Watson, Y. Benjamini, and R.A. Stine. Volume and Area of Oil Fields and Their Impact on Order of Discovery, 1980. http://www.osti.gov/scitech/servlets/purl/6037591.
}

\author{
  Jonathan D. Rosenblatt \email{john.ros.work@gmail.com}
}


\seealso{
  \code{\link{estimate.b.k}}, \code{\link{initializeRdsObject}}, \code{\link{makeRdsSample}}
}

\examples{
# See in estimate.b.k()
}


