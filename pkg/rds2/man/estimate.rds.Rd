\name{estimate.rds}
\alias{estimate.rds}
\alias{estimate.rds.fix.theta}
\title{ ML estimation of population size and degree distribution in respondent driven samples }
\usage{
  estimate.rds(data, Sij, init, const, arc = FALSE,
    maxit = 10000, initial.thetas, theta.minimum,
    theta.range)
    
     estimate.rds.fix.theta(data, Sij, init, const, arc = FALSE,
    maxit = 10000, theta)
}
\arguments{
  \item{data}{An integer vector of the degree of each individual sampled (including zeroes).}

  \item{Sij}{An integer matrix of the counts of individuals sampled with rank i up to sampling period j. Row names are required and represent degrees.}

  \item{init}{Optional initialization values for the parameters except theta (see initial.thetas).}

  \item{const}{Control the concentration around zero of the mapping of the theta parameter from the constrained domain to the real line.}

  \item{arc}{Should the degrees in the snowball be taken in account (TRUE) or just the size of the snowball (FALSE). Following the lines of [see reference] defaults to FALSE.}

  \item{maxit}{Number of maximal optimization iterations.}

  \item{initial.thetas}{ Vector of initial values of theta parameter.}
  
  \item{theta}{ The known value of \eqn{\theta} used by \code{estimate.rds.fix.theta}.}

  \item{theta.minimum}{Lower bound on the allowed values of theta.}

  \item{theta.range}{Allowed range of theta values. From theta.minimum to theta.minimum+range.}

}
\description{
    ML estimation of population size and degree distribution in respondent driven samples. 
}
\value{
A list with the estimates for each starting value of \eqn{\theta}. Each components includes:
	\item{beta}{See in \emph{details}}
	\item{theta}{See in \emph{details}}
	\item{Nj}{See in \emph{details}}
	\item{iterations}{The number of \code{\link{optim}} iterations until convergence.}
	\item{likelihood.optimum}{The (log) likelhood at convergence.}	
}
\details{
\code{estimate.rds} performs maximum likelihood estimation of the population
  size, degree distribution and "coefficient of discoverability" as described in the reference. 
  \code{estimate.rds.fix.theta} performs the same estimation but with \eqn{\theta} assumed known.\cr 
  The method assumes the probability of sampling an individual with degree i at time j is given by: 
  \deqn{\beta I[j] (N[i]-S[i,j]) i^\theta }
  where   \eqn{\beta} is the base rate of sampling, \eqn{\theta} is the "coefficient of discoverability", \eqn{I[j]} it the number of individuals sampled up to time \eqn{j}, 
   \eqn{N[i]} is the number of individuals with degree \eqn{i} in the population, and \eqn{S[i,j]} is the number of individuals with degree \eqn{i} sampled up to time \eqn{j}.
   }
\examples{
data(simulation, package='chords')
temp.data<- unlist(data3[1,7000:7500])
initial.thetas<- c(1,2)
rds.result<- estimate.rds(temp.data, Sij = make.Sij(temp.data), initial.thetas = initial.thetas , 
	arc = FALSE, maxit = 100, const = 0.5, theta.minimum = -0.5, theta.range = 2)
str(rds.result)
plot(rds.result$Nj, type='h', xlab='Degree', ylab=expression(N[j]), main='Estimated Degree Distribution')	
}
\author{
  Jonathan D. Rosenblatt \email{john.ros@gmail.com}
}
\references{
	Respondent Driven Sampling as an Epidemic Process
	Berchenko Y., Rosenblatt J.D., White R.G.,Frost S.D.W. (2012) 
  \emph{####.}  
}


