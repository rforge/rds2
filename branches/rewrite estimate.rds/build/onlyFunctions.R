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
#' @author Jonathan Rosenblatt 
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



#' Make the Sij matrix assuming no droouts from the snowball
#' @param data 
#' @return TBC
#' @author johnros
#' @export
make.Sij<- function(data){
	result<- matrix(0, ncol=length(data), nrow=max(data))
	for (i in 1:(length(data)-1)){ result[data[i],i+1]<- result[data[i],i+1]+1 }  
	result<- t(apply(result, 1 , cumsum))
	non.empty.ind<- apply(result, 1, sum)!=0
	.row.names<- as.character(seq(along.with=non.empty.ind)[non.empty.ind])
	result<- result[non.empty.ind,]
	rownames(result)<- .row.names
	return(result)
}



#' Inverts the map of theta to the real line.  
#' @param qnorm.theta 
#' @param const 
#' @param range 
#' @param minimum 
#' @return TBC
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
#' @return TBC
#' @author johnros
qnorm.theta<- function(theta, const, range=20, minimum=-10){
	normalized.theta<- (theta-minimum)/range
	result<- const*qnorm(normalized.theta)
	return(result)
}

## Test:
# inv.qnorm.theta(10, 100)
# inv.qnorm.theta(10, 10)
# inv.qnorm.theta(20000,10)
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




#' Compares the output of the optimizaton for different initialization values.
#' @param Sij 
#' @return Used for creating the Snowball size given an Sij matrix
#' @author johnros
#' @export
compute.S<- function(Sij){
	result<- apply(Sij, 2, sum)
	non.null.ind<- cumsum(result!=0)
	result[non.null.ind==0]<- 1    
	return(as.integer(result))
}


## TODO: Fix estimation for a single degree!

#' Same as estimate.rds, but for a fixed theta.
#' @param data 
#' @param Sij 
#' @param init 
#' @param const 
#' @param arc 
#' @param maxit 
#' @param theta 
#' @param ... 
#' @return TBC 
#' @author johnros
#' @useDynLib rds2
#' @export 
#' @examples
#' data(brazil)
#' estimate.rds2(data=data.degree, , Sij = data.Sjt, const=50, arc=FALSE, maxit=1000, theta = 1)
estimate.rds2<- function (data, Sij, init, const, arc=FALSE, maxit=10000, theta, ...) {
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



#' Main function in rds2 package. Returns the ML estimate of population size and degree distribution.
#' 
#' Performs maximum likelihood estimation of the population size, degree distribution and theta [explain] allowing (a) "soft" constraints on theta and (b) an input Sij matrix. 
#' 
#' @param data 
#' @param Sij 
#' @param init 
#' @param const 
#' @param arc 
#' @param maxit 
#' @param initial.thetas 
#' @param theta.minimum 
#' @param theta.range 
#' @param ... 
#' @author Jonathan Rosenblatt
#' @useDynLib rds2
#' @export
#' @examples
#' data(simulation)
#' estimate.rds3(temp.data, Sij = make.Sij(temp.data), initial.thetas = c(1,10), arc = FALSE, maxit = 1000, const = 0.5, theta.minimum = -0.5, theta.range = 2)
estimate.rds3<- function (data, Sij, init, const, arc=FALSE, maxit=10000, initial.thetas, theta.minimum, theta.range, ...) {
	# Look for degrees in the data, so their estimates are non vanishing
	N.j<- rep(0, max(data)) 
	uniques<- unique(data)
	uniques<- uniques[uniques!=0]
	sorted.uniques<- sort(uniques)
	param.size<- length(uniques)
	data.table<-table(data)[-1]
	
	# Compute the size of the snowball along the sample:
	S<- compute.S(Sij)
	
	
	# Wrapper to the likelihood function. Implements constraints on parameters by taking real valued parameters and remapping them.
	likelihood.wrap<- function(par){
		final.result<- -Inf # initialize output
		beta<- exp(par[1])
		theta<- inv.qnorm.theta(par[2], const = const)		
		## TODO: A) Add estimatino of theta. 
		N.j<- rep(0, max(data))
		N.j[sorted.uniques]<- exp(tail(par,-2)) # fill non trivial Nj estimates.
		if(any( N.j[sorted.uniques] < data.table )) return(final.result) # checks that given Njs correspond to estimatable values.		
		
		if(is.numeric(beta) && is.numeric(theta) && !is.infinite(theta) && !is.infinite(beta) ) {
			result<- .C("likelihood", 
					sample=as.integer(data), 
					Sij=as.integer(as.matrix(Sij)),
					S=as.integer(S),				  
					c=as.numeric(beta), 
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
	
	## Initialize estimation:
	
	cap<- max(data)
	cap.arcs<-  sum.ranks(data)
	max.cap.arcs<- max(cap.arcs)
	log.beta<- -log(max(sum.ranks(data)*sum(data)) * max(data)) - 0.01 # beta has to be such that all probabilities in the likelihood are between 0 and 1. In particular, for the last subject sampled.	
	
	
	if(missing(init)) {		
		init<- lapply(initial.thetas, function(x) c(log.c=log.beta, qnorm.theta=qnorm.theta(x,const), log(data.table)+1 ))
	}
	
	likelihood.optim<-lapply(init, function(x) {
				try(optim(par=x, fn=likelihood.wrap, control=list(fnscale=-1, maxit=maxit),...)) 
			}  ) 
	
	prepare.result<- function(x){
		result<- simpleError("Optim did not converge") 
		if(is.list(x)){
			N.j[sorted.uniques]<- exp(tail(x$par,-2))
			result<- list( 
					c=exp(x$par[[1]]), 
					theta=inv.qnorm.theta(x$par[[2]], const=const, range=theta.range, minimum=theta.minimum), 
					Nj=N.j,
					iterations=x$counts,
					likelihood.optimum=x$value)			
		}		
		return(result)
	} 
	## TODO: A) Return only numeric objects
	
	return(lapply(likelihood.optim, prepare.result))							
}



