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
#' Assumes input is a table of degrees and degree counts.
#' @param data The sample
#' @return The sum of degrees for each rank 
#' @author Jonathan Rosnblatt
#' @export
ranks.sum<- function(data){
	data.table<- table(data)
	result<- as.numeric(rownames(data.table))*c(data.table)
	return(result)
}



#' Make the Sij matrix assuming no droouts from the snowball
#' @param data 
#' @return TBC
#' @author johnros
#' @export
make.Sij<- function(data){
	stopifnot(length(unique(data))>1L)
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
#' @param theta.range 
#' @param theta.minimum 
#' @return TBC
#' @author johnros
inv.qnorm.theta<- function(qnorm.theta, ...){
#	normalized.theta<- pnorm(qnorm.theta / const)
#	normalized.theta2<- (normalized.theta*theta.range) + theta.minimum
#	return(normalized.theta2)
#	exp(qnorm.theta)
	qnorm.theta
}



#' Maps the theta parameter to the real line.
#' @param theta 
#' @param const Controls the 
#' @param theta.range 
#' @param theta.minimum 
#' @return TBC
#' @author johnros
qnorm.theta<- function(theta, ...){
#	normalized.theta<- (theta-theta.minimum)/theta.range
#	result<- const*qnorm(normalized.theta)
#	return(result)
#	log(theta)
	theta
}





## TODO: B) Function no longer needed as only the best result is returned by estimate.rds

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








generate.rds.control<- function(maxit=2000){
	return(list(maxit=maxit))
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
#'data(simulation, package='chords')
#'temp.data<- unlist(data3[1,7000:7500])
#'(rds.result<- estimate.rds(sampled.degree.vector=temp.data , Sij=make.Sij(temp.data), initial.values=list(theta=c(1,2)), control=generate.rds.control(maxit=20)))
#'plot(rds.result$Nj, type='h', xlab='Degree', ylab=expression(N[j]), main='Estimated Degree Distribution')	
#'var.theta(temp.data, Njs=rds.result$Nj, theta=rds.result$theta, beta=rds.result$beta)

estimate.rds<- function (sampled.degree.vector, Sij, method="BFGS", initial.values, arc=FALSE, control=generate.rds.control(), all.solutions=FALSE) {  	
  	# Initializing:
	max.observed.degree<- max(sampled.degree.vector)
	N.j<- rep(0, max.observed.degree)
	final.result<- NA
	# Look for the degrees in the datawith non trivial estimates:
	Observed.Njs<-table(sampled.degree.vector)[-1] # table of degree counts withuot zero
	Observed.js<- as.numeric(names(Observed.Njs))
	maximal.degree.count<- max(Observed.Njs)
	the.call<- sys.call()
	
	
	# Compute the size of the snowball along the sample:
	S<- compute.S(Sij)
	
	
	# Wrap the likelihood function. Implements constraints on parameters by taking real valued parameters and remapping them to the constrained space.
	likelihood.wrap<- function(par){
		# Initialize:
		likelihood.result<- -9999999 
			
		beta<- exp(par[['canonical.beta']]) # assumes beta is non negative and given in log scale
		theta<- do.call(inv.qnorm.theta, c(qnorm.theta=par[['canonical.theta']], control))		
		N.j<- rep(0, max.observed.degree)
		## Does this sbsetting work when initial values are not specified? 
		observed.js.indexes<- 	sapply(Observed.js, function(x) grep(paste("logNjs.",x,"$",sep=""), names(par) )   )
		N.j[Observed.js]<- exp(par[observed.js.indexes]) # fill non trivial Nj estimates.
		
		# Checking estimates are within allowed range:
		if(	isTRUE( any( round(N.j[Observed.js],10) < round(Observed.Njs,10) ) || beta >= 1 / ( sum(N.j) * max(N.j) * max.observed.degree^theta ))){
			return(likelihood.result)
		} 
		
		if(is.numeric(beta) && is.numeric(theta) && !is.infinite(theta) && !is.infinite(beta) ) {
			result<- .C("likelihood", 
					sample=as.integer(sampled.degree.vector), 
					Sij=as.integer(as.matrix(Sij)),
					S=as.integer(S),				  
					beta=as.numeric(beta), 
					theta=as.numeric(theta), 
					Nj=as.numeric(N.j), 
					observed_degrees=as.integer(rownames(Sij)),
					n=as.integer(length(sampled.degree.vector)), 
					N=as.integer(length(N.j)),
					N_observed=as.integer(nrow(Sij)),				  
					arc=arc,
					result=as.double(0))		  
			
		}	
		if(!is.infinite(result$result)) likelihood.result<-result$result
		return(likelihood.result)
		
	}
	
	
	
	
	
	# TODO: B) Automate initialization values for theta. (maybe using moments of intervals between samples).
	# TODO: A) Use better N_j initialization: IDW or inflated values.
	generate.initial.values<- function(initial.values){
		# compute maximal beta given other parameter values:
		initial.log.beta.function<- function(theta) {
			beta<- 1/ (sum(sampled.degree.vector) * maximal.degree.count * max.observed.degree^theta )
			return(log(beta) - 0.01)
		}
		
		
		
			
		# initiate with only theta given:
		if(length(initial.values)==1L){
			stopifnot(is.numeric(initial.values$theta))			
			wrap.initial.values<- function(theta) {
				c(
						canonical.beta=initial.log.beta.function(theta), 
						canonical.theta=do.call(qnorm.theta, c(theta=theta, control)),
						logNjs= log(Observed.Njs)+0.1 
						)
			}
			returned.initial.values<- lapply(initial.values$theta, wrap.initial.values)	
		}
		
		
		
		# initiate with all parameter given:
		else if(length(initial.values)==3L){
			stopifnot(all(is.numeric(initial.values$theta), is.numeric(initial.values$beta), sapply(initial.values$Njs, is.numeric)))
			wrap.initial.values<- function(beta, theta, Njs){
				c(
						canonical.beta=log(beta) ,
						canonical.theta=do.call(qnorm.theta, c(theta=theta, control)) ,
						logNjs=log(Njs)#+1
				)
			}
			returned.initial.values<- mapply(wrap.initial.values, beta= initial.values$beta, theta=initial.values$theta, Njs=initial.values$Njs, SIMPLIFY=FALSE )
		}
		
		
		# initiate with theta and beta given:
		else if(length(initial.values)==2L){
			stopifnot(all(is.numeric(initial.values$theta), is.numeric(initial.values$beta)))
			wrap.initial.values<- function(beta, theta){				
				output<- c(
						canonical.beta=log(beta) ,
						canonical.theta=qnorm.theta(theta) ,
						logNjs=log(Observed.Njs)+0.1
				)				
				log.beta.limit<- initial.log.beta.function(theta)
				
				if(beta > exp(log.beta.limit)){
					message(paste('Initial beta too large. Forcing beta=',signif(exp(log.beta.limit),4), " instead of beta=",signif(beta,4), sep=""))
					output[['canonical.beta']]<- log.beta.limit
				}				
				return(output)
			}
			
			# take a the list of beta and thetas and return a list with initialization values
			returned.initial.values<-list()
			for (theta in initial.values$theta){
				for (beta in initial.values$beta){
					returned.initial.values<- c(returned.initial.values, list(wrap.initial.values(beta, theta)))
				}
			}				
		}	
		
		
		return(returned.initial.values)
	}
		
	
	
	
	
	
	
	## Initialize estimation:	 	
	initial.values.list<- generate.initial.values(initial.values)				
		
	# In case of convergence problems, try different optimization methods. In particular: "Nelder-Mead", "SANN", 
	optim.wrap<- function(x) try(optim(par=x, fn=likelihood.wrap, method=method , control=list(fnscale=-1, maxit=control$maxit)))
	
	likelihood.optim<-lapply(initial.values.list,  optim.wrap )
	
	
	prepare.result<- function(x){
		result<- simpleError("Optim did not converge") 
		if(length(x)==5L){
			observed.js.indexes<- sapply(Observed.js, function(y) grep(paste("logNjs.",y,"$",sep=""), names(x$par)))
			N.j[Observed.js]<- exp(x$par[observed.js.indexes]) # fill non trivial Nj estimates.		
			
			result<- list( 
					beta=exp(x$par[['canonical.beta']]), 
					theta=do.call(inv.qnorm.theta, c(qnorm.theta=x$par[['canonical.theta']], control)),
					initial.values=initial.values,
					Nj=N.j,
					iterations=x$counts,
					likelihood.optimum=x$value, 
					call=the.call)			
		}		
		return(result)
	} 

	temp.result<- lapply(likelihood.optim, prepare.result)
	
	 if(  any(sapply(temp.result, length) > 2 )   ){
			clean.temp.result<- as.list(temp.result[sapply(temp.result, length) > 2])
			if(!all.solutions){ 
				final.result<- temp.result[[  which.max(sapply(clean.temp.result, function(x) x$likelihood.optimum)) ]]
			}
			else{
				final.result<- clean.temp.result
			} 
		}
		else{
			message('Estimation did not converge. Try differet intialization values or optimization method.')
		}
		
		return(final.result)							
}

##### Testing: 
#data(simulation, package='chords')
#temp.data<- unlist(data3[1,7000:7500])
## Initialize only with thetas:
#(rds.result<- estimate.rds(sampled.degree.vector=temp.data , Sij=make.Sij(temp.data), method="BFGS", initial.values=list(theta=c(1,2,10)), 
#					control=generate.rds.control()))
#str(rds.result)
#plot(rds.result$Nj, type='h', xlab='Degree', ylab=expression(N[j]), main='Estimated Degree Distribution')	
#var.theta(temp.data, Sij = make.Sij(temp.data), Njs=rds.result$Nj, theta=rds.result$theta, 
#		beta=rds.result$c)
#
#			 
#
#createDegreeCount <- function(network.size, network.density) {
#	neighbour.count<- lapply(seq(1, network.size), function(x) rbinom(1, network.size, network.density))
#	return(table(unlist(neighbour.count[neighbour.count>0])))
#}
## Good values1: network.size = 500, network.density = 0.1, theta<- 3, beta<- 1e-10
##(Njs<- createDegreeCount(network.size = 1e5, network.density = 0.01))
## Manual creation:
#Njs<- c(100,100,100,100); names(Njs)<- c("1","50","100","1000"); Njs<- as.table(Njs)
## Good values2: theta<- 3, beta<- 5e-10
#theta<- 1
#beta<- 5e-8
#tail(degree.sampled.vec<- generate.sample(theta, Njs, beta, sample.length=1e3))
#plot(degree.sampled.vec, type='h', main='Sample')
#
#
## Testing with initialization from **true** values:
#str(rds.result<- estimate.rds(degree.sampled.vec, Sij = make.Sij(degree.sampled.vec), method='BFGS',				
#				initial.values = list(theta=c(theta), beta=c(beta), Njs=list(Njs)),  
#				control = generate.rds.control( maxit = 3000)))
#str(rds.result)
#plot(Njs, type='h', lwd=2.5)	
#points(rds.result$Nj, type='h', col='orange', lwd=2)
#points(rds.result[[1]]$Nj, type='h', col='orange', lwd=2)
#points(rds.result[[2]]$Nj, type='h', col='orange', lwd=2)
#points(rds.result[[3]]$Nj, type='h', col='orange', lwd=2)
#
#
## Initializing with theta only. (true values unknown):
#str(rds.result<- estimate.rds(degree.sampled.vec, Sij = make.Sij(degree.sampled.vec), arc=FALSE, 
#				initial.values = list(theta=c(-1, 0.5, 2, 8)), method="BFGS",
#				control = generate.rds.control(maxit = 1000)))
#str(rds.result)
#plot(Njs, type='h', lwd=2.5); points(rds.result[[1]]$Nj, type='h', col='orange', lwd=2)
#plot(Njs, type='h', lwd=2.5); points(rds.result[[2]]$Nj, type='h', col='orange', lwd=2)
#plot(Njs, type='h', lwd=2.5); points(rds.result[[3]]$Nj, type='h', col='orange', lwd=2)
#plot(Njs, type='h', lwd=2.5); points(rds.result[[4]]$Nj, type='h', col='orange', lwd=2)
#
#
## Initializing with thetas and betas (unknown values):
#str(rds.result<- estimate.rds(degree.sampled.vec, Sij = make.Sij(degree.sampled.vec), arc=FALSE, 
#				initial.values = list(theta=c(-1, 0.5, 2, 8), beta=c(1e-13, 1e-1, 1e-5)), 
#				method="BFGS", control = generate.rds.control(maxit = 1000)))
#str(rds.result)









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





var.theta<- function(sampled.degree.vector, Sij=make.Sij(sampled.degree.vector), Njs, theta, beta){
	I<- compute.S(Sij)
	N<- sum(Njs*seq(along.with=Njs))
	Nj.uniques<- sort(unique(sampled.degree.vector))
	total.sampling<- length(sampled.degree.vector)
	make.value<- Vectorize(function(k,t){
				if(!(k %in% rownames(Sij))) return(0)
				else return(
							log(k)^2* k^theta* I[t]/N *(Njs[k]/N - Sij[as.character(k),t]/N)
					)
			})
	information<- beta * N * sum(outer(Nj.uniques, seq(1,total.sampling), FUN="make.value"))
	return( (1/information)/N )
}
## Example:
#var.theta(degree.sampled.vec, make.Sij(degree.sampled.vec), rds.result$Nj, rds.result$theta, rds.result$beta )
