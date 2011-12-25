#----------------- Utility functions---------------#

# Old versions #
#load('/home/johnros/Projects/Yakir/oldFunctions.Rdata')
####################################################


#' Estimates oppulation size and degree distribution in respondant driven samples(RDS).
#' 
#' \tabular{ll}{
#' Package: \tab rds2\cr
#' Type: \tab Package\cr
#' Version: \tab 0.5\cr
#' Date: \tab 2011-12-25\cr
#' License: \tab GPL (>= 2)\cr
#' }
#' 
#' Maximum likelihood estimation of population size using the methodology in [cite].
#'
#' @name rds2-package
#' @aliases rds2
#' @docType package
#' @title Population size estimation for despondant driven sampling
#' @author Jonathan Rosenblatt \email{john.ros@gmail.com}
#' @references
#' @keywords 
#' @seealso 
#' @examples
#' 
NA


#' Utility functions for creating snowball matrix.
#' @param data The sample
#' @return The sum of degrees for each rank 
#' @author Jonathan Rosnblatt
#' @export
sum.ranks<- function(data){
  data.table<- table(data)[-1]
  result<- as.numeric(rownames(data.table))*c(data.table)
  return(result)
}

#' Inverts the map of theta to the real line.  
#' @param qnorm.theta 
#' @param const 
#' @param range 
#' @param minimum 
#' @return 
#' @author johnros
inv.qnorm.theta<- function(qnorm.theta, const, range=20, minimum=-10){
	normalized.theta<- pnorm(qnorm.theta / const)
	normalized.theta2<- (normalized.theta*range) + minimum
	return(normalized.theta2)
}



#' Maps the theta parameter to the real line.
#' @param theta 
#' @param const Controls the 
#' @param range 
#' @param minimum 
#' @returnType 
#' @return 
#' @author johnros
qnorm.theta<- function(theta, const, range=20, minimum=-10){
	normalized.theta<- (theta-minimum)/range
	result<- const*qnorm(normalized.theta)
	return(result)
}

## Test:
# inv.qnorm.theta(10, 100)
# inv.qnorm.theta(10, 10)
#inv.qnorm.theta(20000,10)
# inv.qnorm.theta(25,200)
# inv.qnorm.theta(25,200, range=2)
# inv.qnorm.theta(25,200, range=1)
# curve(inv.qnorm.theta(x,const=0.6, range=2, minimum=1), -10,10)
# qnorm.theta(25, range=30, minimum=-1, const=0.1)
# qnorm.theta(theta=-0.5, const=10, range=2, minimum=-1)
# inv.qnorm.theta(qnorm.theta(5,100),100)
# curve(qnorm.theta(x,const=0.6, range=1, minimum=0), -2,2)
# curve(qnorm.theta(x,const=30, range=1, minimum=-1), -2,2)




#' Compares the output of the optimizaton for different initialization values.
#' @param estimation An rds2 estimation object
#' @return Used only for printing
#' @author johnros
#' @export
comparison <- function (estimation) {
  non.na<- !is.na(estimation)
  print(signif(cbind(sapply(estimation[non.na], function(x) return(x$Nj)))))
  print(sapply(estimation[non.na], function(x) sum(x$Nj)))
  print(sapply(estimation[non.na], function(x) x$theta))
  print(sapply(estimation[non.na], function(x) x[[4]]$value))
}









# Optimization stage:

## TODO: Fix estimation for a single degree!

# Note: make sure likelihood.cpp is compiled for the right system!
# For 32 bit machines: dyn.load("/home/johnros/workspace/rds2/pkg/likelihoodVer9.32.so")
# dyn.load("/home/johnros/workspace/rds2/pkg/likelihoodVer9.64.so")





#' Main function in rds2 package. Returns the ML estimate of population size and degree distribution. 
#' @param data 
#' @param init 
#' @param const 
#' @param arc 
#' @param maxit 
#' @param ... 
#' @returnType 
#' @return 
#' @author Jonathan Rosenblatt
#' @export
estimate.rds<- function (data, init, const, arc, maxit=20000, ...) {
  # Look for degrees in the data, so their estimates are nony vanishing
  N.j<- rep(0, max(data)) 
  uniques<- unique(data)
  uniques<- uniques[uniques!=0]
  s.uniques<- sort(uniques)
  param.size<- length(uniques)
  data.table<-table(data)[-1]
  
  likelihood.wrap1<- function(par){
    final.result<- -Inf
	  ccc<- exp(par[1])
	  theta<- inv.qnorm.theta(par[2], const = const)
	  N.j<- rep(0, max(data))
	  N.j[s.uniques]<- exp(tail(par,-2))
	  if(any( N.j[s.uniques] < data.table )) return(final.result)	  
	  
	  if(is.numeric(ccc) && is.numeric(theta) && !is.infinite(theta) && !is.infinite(ccc) ) {
		  result<- .C("likelihood", 
				  sample=as.integer(data), 
				  c=as.numeric(ccc), 
				  theta=as.numeric(theta), 
				  Nj=as.numeric(N.j), 
				  constant=as.numeric(const),
				  n=as.integer(length(data)), 
				  N=as.integer(length(N.j)),
				  arc=arc,
				  result=as.double(0))
		  final.result<- result$result
	  }	  
	  #print(c(as.integer(data),log.c=as.numeric(par[1]),qnorm.theta=as.numeric(par[2]),Nj=as.numeric(N.j),constant=as.numeric(const)))
	  
	  return(final.result)
  }
  cap<- max(data)
  cap.arcs<-  sum.ranks(data)
  max.cap.arcs<- max(cap.arcs)
  log.c<- -log(max(sum.ranks(data)*sum(data)) * max(data)) -1
  
  if(missing(init)) { 
	  init<- list( 
			  six0.2= c(log.c=log.c, qnorm.theta=qnorm.theta(0.2,const), log(data.table)+1 ),
			  #six0.5= c(log.c=log.c, qnorm.theta=qnorm.theta(0.5,const), log(data.table)+1 ),
			  #six0.9= c(log.c=log.c, qnorm.theta=qnorm.theta(0.9,const), log(data.table)+1 ),
			  six1= c(log.c=log.c, qnorm.theta=qnorm.theta(1,const), log(data.table)+1 ),
			  #six1.1= c(log.c=log.c, qnorm.theta=qnorm.theta(1.1,const), log(data.table)+1 ),
			  #six1.5= c(log.c=log.c, qnorm.theta=qnorm.theta(1.5,const), log(data.table)+1 ),
			  six2= c(log.c=log.c, qnorm.theta=qnorm.theta(2,const), log(data.table)+1 )			  
	  )}
  
  likelihood.optim<-lapply(init, function(x) {
			  try(optim(par=x, fn=likelihood.wrap1, control=list(fnscale=-1, maxit=maxit),...)) }  )  
  
  prepare.result<- function(x){    
    N.j[s.uniques]<- exp(tail(x$par,-2))
    result<- list( 
      c=exp(x$par[[1]]), 
      theta=inv.qnorm.theta(x$par[[2]],const), 
      Nj=N.j, 
	  x  )
    return(result)
  }  
  return(lapply(likelihood.optim, function(x) try(prepare.result(x))))
}




#' Main function in rds2 package. Returns the ML estimate of population size and degree distribution.
#' @param data 
#' @param Sij 
#' @param init 
#' @param const 
#' @param arc 
#' @param maxit 
#' @param theta 
#' @param ... 
#' @returnType 
#' @return 
#' @author johnros
#' @export
#' @example 
#' data(Cornell)
#' Cornell<- estimate.rds2(data=data.degree, , Sij = data.Sjt, const=50, arc=FALSE, maxit=10000)
#' comparison(Cornell)
estimate.rds2<- function (data, Sij, init, const, arc, maxit=10000, theta, ...) {
  # Look for degrees in the data, so their estimates are nony vanishing
  N.j<- rep(0, max(data)) 
  uniques<- unique(data)
  uniques<- uniques[uniques!=0]
  s.uniques<- sort(uniques)
  param.size<- length(uniques)
  data.table<-table(data)[-1]
  
  S<- compute.S(Sij)
  
  likelihood.wrap1<- function(par){
    final.result<- -Inf
	  ccc<- exp(par[1])
	  N.j<- rep(0, max(data))
	  N.j[s.uniques]<- exp(tail(par,-1))
	  if(any( N.j[s.uniques] < data.table )) return(final.result)	  
	  
	  if(is.numeric(ccc) && is.numeric(theta) && !is.infinite(theta) && !is.infinite(ccc) ) {
		  result<- .C("likelihood", 
				  sample=as.integer(data), 
				  Sij=as.integer(as.matrix(Sij)),
				  S=as.integer(S),				  
				  c=as.numeric(ccc), 
				  theta=as.numeric(theta), 
				  Nj=as.numeric(N.j), 
				  constant=as.numeric(const),
				  observed_degrees=as.integer(rownames(Sij)),
				  n=as.integer(length(data)), 
				  N=as.integer(length(N.j)),
				  N_observed=as.integer(nrow(Sij)),				  
				  arc=arc,
				  result=as.double(0))		  
		  final.result<- result$result
	  }	  	  
	  return(final.result)
  }
  cap<- max(data)
  cap.arcs<-  sum.ranks(data)
  max.cap.arcs<- max(cap.arcs)
  log.c<- -log(max(sum.ranks(data)*sum(data)) * max(data)) -1
  
  if(missing(init)) { 
	  init<- list( 
			  six= c(log.c=log.c, log(data.table)+1 )

	  )}
  
   likelihood.optim<-lapply(init, function(x) {
			  try(optim(par=x, fn=likelihood.wrap1, control=list(fnscale=-1, maxit=maxit),...)) 
		  }  ) 
  
  prepare.result<- function(x){    
    N.j[s.uniques]<- exp(tail(x$par,-1))
    result<- list( 
      c=exp(x$par[[1]]), 
      theta=theta,
      Nj=N.j, 
	  x  )
    return(result)
  }  
  return(lapply(likelihood.optim, function(x) try(prepare.result(x))))
}


