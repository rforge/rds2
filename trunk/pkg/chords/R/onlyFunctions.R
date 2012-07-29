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
#' @return TBC 
#' @author johnros
#' @useDynLib rds2
#' @export 
#' @examples
#' data(brazil)
#' estimate.rds2(data=data.degree, , Sij = data.Sjt, const=50, arc=FALSE, maxit=1000, theta = 1)
estimate.rds.fix.theta<- function (data, Sij, init, const, arc=FALSE, maxit=10000, theta) {
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
		beta<- exp(par[1])
		N.j<- rep(0, max(data))
		N.j[s.uniques]<- exp(tail(par,-1))
		if(any( N.j[s.uniques] < data.table )) return(final.result)	  
		
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
	cap<- max(data)
	cap.arcs<-  sum.ranks(data)
	max.cap.arcs<- max(cap.arcs)
	log.beta<- -log(max(sum.ranks(data)*sum(data)) * max(data)) -1
	
	if(missing(init)) { 
		init<- list( 
				six= c(log.c=log.beta, log(data.table)+1 )
		
		)}
	
	likelihood.optim<-lapply(init, function(x) {
				try(optim(par=x, fn=likelihood.wrap1, control=list(fnscale=-1, maxit=maxit))) 
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
#' @author Jonathan Rosenblatt
#' @useDynLib rds2
#' @export
#' @examples
#' data(simulation)
#' estimate.rds3(data= temp.data, Sij = make.Sij(temp.data), initial.thetas = c(1,10), arc = FALSE, maxit = 1000, const = 0.5, theta.minimum = -0.5, theta.range = 2)

estimate.rds<- function (data, Sij, init, const=0.5, arc=FALSE, maxit=10000, initial.thetas=c(0.5, 1, 1.5), theta.minimum=-1, theta.range=2) {
	# Verifying input:
  stopifnot(!all(missing(x=initial.thetas), missing(x=theta.minimum), missing(x=theta.range)))
  
  	
  	# Initializing:
	# Look for degrees in the data, so their estimates are non vanishing
	N.j<- rep(0, max(data)) 
	uniques<- unique(data)
	uniques<- uniques[uniques!=0]
	sorted.uniques<- sort(uniques)
	param.size<- length(uniques)
	data.table<-table(data)[-1]
	final.result<- NA
	
	# Compute the size of the snowball along the sample:
	S<- compute.S(Sij)
	
	
	# Wrapper to the likelihood function. Implements constraints on parameters by taking real valued parameters and remapping them.
	likelihood.wrap<- function(par){
		final.result<- -Inf # initialize output
		beta<- exp(par[1])
		theta<- inv.qnorm.theta(par[[2]], const=const, range=theta.range, minimum=theta.minimum)		
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
				try(optim(par=x, fn=likelihood.wrap, control=list(fnscale=-1, maxit=maxit))) 
			}  ) 
	
	prepare.result<- function(x){
		result<- simpleError("Optim did not converge") 
		if(is.list(x)){
			N.j[sorted.uniques]<- exp(tail(x$par,-2))
			result<- list( 
					beta=exp(x$par[[1]]), 
					theta=inv.qnorm.theta(x$par[[2]], const=const, range=theta.range, minimum=theta.minimum), 
					Nj=N.j,
					iterations=x$counts,
					likelihood.optimum=x$value)			
		}		
		return(result)
	} 

	temp.result<- lapply(likelihood.optim, prepare.result)
	
	if(any(sapply(temp.result, length)==5L)) {
		clean.temp.result<- temp.result[sapply(temp.result, length)==5L]
		final.result<- temp.result[[which.max(sapply(clean.temp.result, function(x) x$likelihood.optimum))]]
	}
	else{
		message('Estimation did not converge. Try differet intialization values')
	}
	
	return(final.result)							
}

data(simulation, package='chords')
temp.data<- unlist(data3[1,7000:7500])
(rds.result<- estimate.rds(temp.data, Sij = make.Sij(temp.data), initial.thetas = c(0.5,1,2, 100) , 
		arc = FALSE, maxit = 100, const = 0.5, theta.minimum = -0.5, theta.range = 2))










generate.sample <- function(theta, Njs, beta, sample.length, double.sampling.warning.thresh=0.05) {
	# Initializing:
	js<- as.numeric(names(Njs))
	N<- sum(Njs)
	maximal.beta<- 1/ (N * max(Njs) * max(js)^theta )
	double.sampling.warnings<- FALSE
	
	
	# Input validations:
	error.messge<- paste("beta value too large. Errors might occur. \nFor the given degree frequencies, consider keeping it under ", signif(maximal.beta), sep="" )
	if(beta > maximal.beta) message(error.messge) # Note: larger beta favour sampling non 0 degrees.
	stopifnot(min(js)>0)
	
	
	
	# Preparing to sample:
	Uik<- rep(0, length.out=max(js))
	names(Uik)<- seq(along.with=Uik)
	snowball<- 1
	degree.sampled.vec<- rep(NA, sample.length)
	# Start sampling:
	for(i in seq(1, sample.length)){	 
		pi.function<- function(k) beta * k^theta * snowball * (Njs[paste(k)] - Uik[k] )
		n.pi.function<- function(k) 1-pi.function(k)
		# Compute the probabilities of sampling degree k=0:
		no.sample.probability<- prod(unlist(lapply(js, n.pi.function )))
		sample.someone<- as.logical(rbinom(1, 1, prob = 1 - no.sample.probability ))		
		if(sample.someone){
			# Compute the probabilities of sampling degree k>0:
			probs.function<- function(k){
				nominator<-  pi.function(k) * prod(unlist(lapply(js[js!=k], n.pi.function)))  
				# Note: the true probability is not needed due to the normalization in sample(). A constant is used insted.
				# denominator<-  1 - prod(unlist(lapply(js, n.pi.function)))
				denominator<-  1 			
				return(nominator/denominator)		
			}
			degree.probabilities<- unlist(lapply(js, probs.function))
			if( sum(degree.probabilities)+ no.sample.probability <  1- double.sampling.warning.thresh) double.sampling.warnings<- TRUE				
			degree.sampled<- sample(x=js, prob=degree.probabilities, size=1)      
			Uik[degree.sampled]<- Uik[degree.sampled]+1
			snowball<- sum(Uik)
			degree.sampled.vec[i]<- degree.sampled    
		}  
		else{ 
			degree.sampled.vec[i]<- 0
		}
	}
	
	# Exiting:
	if(double.sampling.warnings) message("Double sampling probabilities non-negligeable. Consider lower beta values.")
	return(degree.sampled.vec)
}





var.theta<- function(data, Sij, Njs, theta, beta){
	I<- compute.S(Sij)
	N<- sum(Njs*seq(along.with=Njs))
	Nj.uniques<- sort(unique(data))
	total.sampling<- length(data)
	make.value<- Vectorize(function(k,t){
				if(!(k %in% rownames(Sij))) return(0)
				else return(
							log(k)^2* k^theta* I[t]/N *(Njs[k]/N - Sij[as.character(k),t]/N)
					)
			})
	information<- beta * N * sum(outer(Nj.uniques, seq(1,total.sampling), FUN="make.value"))
	return( (1/information)/N )
}
