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



# Transforms beta back and forth to canonical (real) scale:
generate.beta.scale<- function(scale=10e10) return(scale)


transform.beta <- function(beta, scale=generate.beta.scale()){
	canonical.beta<- log(beta*scale)
	return(canonical.beta)
}




inv.transform.beta<- function(canonical.beta, scale=generate.beta.scale()){
	beta<- exp(canonical.beta)/scale
	return(beta)
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








generate.rds.control<- function(
		maxit=500, 
		method="Nelder-Mead",
		Nj.inflations=c(1,2,10),
		beta.inflations= exp(-seq(2,0, by=-0.02)),
		fnscale=1000,
		){
	return(list(maxit=maxit, method=method, Nj.inflations=Nj.inflations, beta.inflations=beta.inflations, fnscale=fnscale))
}
## Testing:
#generate.rds.control()




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

estimate.rds<- function (sampled.degree.vector, Sij, initial.values, arc=FALSE, control=generate.rds.control(), all.solutions=FALSE) {  	
  	# Initializing:
	max.observed.degree<- max(sampled.degree.vector)
	N.j<- rep(0, max.observed.degree)
	final.result<- NA
	# Look for the degrees in the data with non trivial estimates:
	Observed.Njs<-table(sampled.degree.vector)[-1] # table of degree counts withuot zero
	Observed.js<- as.numeric(names(Observed.Njs))
	maximal.degree.count<- max(Observed.Njs)
	the.call<- sys.call()
	method <- control$method
	maxit <- control$maxit
	fnscale <- control$fnscale
	
	
	# Compute the size of the snowball along the sample:
	S<- compute.S(Sij)
	
	
	# Wrap the likelihood function. Implements constraints on parameters by taking real valued parameters and remapping them to the constrained space.
	likelihood.wrap<- function(par){
		# Initialize:
		likelihood.result<- -9999999 
			
		beta<- inv.transform.beta(par[['canonical.beta']]) # assumes beta is non negative and given in log scale
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
			return(transform.beta(beta))
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
						canonical.beta=transform.beta(beta) ,
						canonical.theta=do.call(qnorm.theta, c(theta=theta, control)) ,
						logNjs=log(Njs)
				)
			}
			returned.initial.values<- mapply(wrap.initial.values, beta= initial.values$beta, theta=initial.values$theta, Njs=initial.values$Njs, SIMPLIFY=FALSE )
		}
				
		# initiate with theta and beta given:
		else if(length(initial.values)==2L){
			stopifnot(all(is.numeric(initial.values$theta), is.numeric(initial.values$beta)))
			wrap.initial.values<- function(beta, theta){				
				output<- c(
						canonical.beta=transform.beta(beta) ,
						canonical.theta=qnorm.theta(theta) ,
						logNjs=log(Observed.Njs)+0.1
				)				
				log.beta.limit<- initial.log.beta.function(theta)
				
				if(beta > inv.transform.beta(log.beta.limit)){
					message(paste('Initial beta too large. Forcing beta=',signif(inv.transform.beta(log.beta.limit),4), " instead of beta=",signif(beta,4), sep=""))
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
	optim.wrap<- function(x) try(optim(par=x, fn=likelihood.wrap, method=method , control=list(fnscale=-1, maxit=maxit, fnscale=fnscale)))
	
	likelihood.optim<-lapply(initial.values.list,  optim.wrap )
	
	
	prepare.result<- function(x){
		result<- simpleError("Optim did not converge")
		
		if(length(x)==5L){
			observed.js.indexes<- sapply(Observed.js, function(y) grep(paste("logNjs.",y,"$",sep=""), names(x$par)))
			N.j[Observed.js]<- exp(x$par[observed.js.indexes]) # fill non trivial Nj estimates.		
			
			result<- list( 
					beta=inv.transform.beta(x$par[['canonical.beta']]), 
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
#require(chords)
#data(simulation, package='chords')
#temp.data<- unlist(data3[1,7000:7500])
## Initialize only with thetas:
#(rds.result<- estimate.rds(sampled.degree.vector=temp.data , Sij=make.Sij(temp.data), initial.values=list(theta=c(1,2,10)), 
#					control=generate.rds.control(maxit=4, method="BFGS")))
#str(rds.result)
#plot(rds.result$Nj, type='h', xlab='Degree', ylab=expression(N[j]), main='Estimated Degree Distribution')	
#var.theta(temp.data, Sij = make.Sij(temp.data), Njs=rds.result$Nj, theta=rds.result$theta, 
#		beta=rds.result$c)




			 

createDegreeCount <- function(network.size, network.density) {
	neighbour.count<- lapply(seq(1, network.size), function(x) rbinom(1, network.size, network.density))
	return(table(unlist(neighbour.count[neighbour.count>0])))
}
# Good values1: network.size = 500, network.density = 0.1, theta<- 3, beta<- 1e-10
#(Njs<- createDegreeCount(network.size = 1e5, network.density = 0.01))
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








#### Grid Search ####


# Creates loose bound on maximal (log) beta values given theta and data:
initialLogBetas <- function(sampled.degree.vector, theta) {
	max.observed.degree<- max(sampled.degree.vector)
	maximal.degree.count<- max(table(sampled.degree.vector))
	beta<- 1/ (sum(sampled.degree.vector) * maximal.degree.count * max.observed.degree^theta )
	return(log(beta)+log(0.9))
}
## Testing:
#initialLogBetas(sampled.degree.vector = rpois(200,2), theta = 1)





# Given a list of thetas and data, guesses beta values for each theta:
makeBetaGrid<- function(sampled.degree.vector, thetas, betas, beta.inflations){
	result<-list()
	for (theta in thetas){
		for(beta in betas){
			for (inflation in beta.inflations){
				#max.beta<- exp(initialLogBetas(sampled.degree.vector, theta))
				result<- c(result, list(c(theta=theta, beta=beta*inflation)))			
			}		
		}
	}
	return(result)
}
## Testing:
#detach("package:chords")
#require(chords)
#Njs<- c(100,100,100,100); names(Njs)<- c("1","50","100","1000"); Njs<- as.table(Njs)
#theta<- 1
#beta<- 2.5e-8
#tail(degree.sampled.vec<- generate.sample(theta, Njs, beta, sample.length=1e4))
#makeBetaGrid(degree.sampled.vec, theta+c(-0.5, 0, 0.5), deflations= c(1,10,100))



# Make grid of Njs given theta and beta:
makeNjGrid<- function(sampled.degree.vector, theta.beta.grid, Nj.inflations){
	# Empirical counts
	Observed.Njs <- table(sampled.degree.vector)[-1] # table of degree counts withuot zero
	Observed.js <- as.numeric(names(Observed.Njs))
	
	# Inflate empirical counts (arbitrary)
	result<- list()
	for (theta.beta in theta.beta.grid){
		for (inflation in Nj.inflations){
			Njs<-  ceiling(Observed.Njs*inflation)
			names(Njs)<- Observed.js
			Njs.list<- list(Njs)
			result<- c(result, list(c(theta.beta, Nj=Njs.list)))
		}
	}
			
	# Return list of options:
	return(result)
}
## Testing:
#makeNjGrid(degree.sampled.vec, theta.beta.grid, c(1,1000))





## Grid search given all three parameters:
grid.all.parameters<- function(sampled.degree.vector, Sij, grid.values, control){
	# Initializing:
	
	temp.estimate.rds<- function(parameters) {
		result<-try(estimate.rds(
						sampled.degree.vector, 
						Sij = Sij,
						control = control,
						initial.values = list(theta=c(parameters$theta), 
								beta=c(parameters$beta), 
								Njs=list(parameters$Nj))), 
				silent=TRUE)
		if(length(result)>2) return(result$likelihood.optimum)
		else(return(-999999))
	}
	likelihoods<-lapply(grid.values, temp.estimate.rds)
		
	# Return most likely solutions
	max.index<- which.max(unlist(likelihoods))
	result<- list(
			theta=grid.values[[max.index]]$theta, 
			beta=grid.values[[max.index]]$beta, 
			Njs=grid.values[[max.index]]$Nj, 
			likelihood=likelihoods[[max.index]]
	) 
	
	return(result)
}
## Testing:
# Create data: 
#detach("package:chords")
#require(chords)
#Njs<- c(20,20,20,20); names(Njs)<- c("1","50","100","1000"); Njs<- as.table(Njs)
#theta<- 1.2
#beta<- 1.2e-8
#tail(degree.sampled.vec<- generate.sample(theta, Njs, beta, sample.length=1e4))
#plot(degree.sampled.vec, type='h', main='Sample')
## Veryfying likelihood computation:
#beta.deflations<- c(1)
#Nj.infaltions<- c(1,2,10)
#thetas<- theta * c(0.9,1,1.1)
#betas<- beta *c(0.9,1,1.1)
#theta.beta.grid<- makeBetaGrid(degree.sampled.vec, thetas, betas, beta.deflations)
#grid.values<- makeNjGrid(degree.sampled.vec, theta.beta.grid, Nj.infaltions)		
#
#grid.all.parameters(degree.sampled.vec, Sij = make.Sij(degree.sampled.vec),				
#		grid.values = grid.values,  
#		control = generate.rds.control( maxit = 20, method="Nelder-Mead"))
#
#estimate.rds(degree.sampled.vec, Sij = make.Sij(degree.sampled.vec), control = generate.rds.control(maxit=10, method="Nelder-Mead"),
#		initial.values = list(theta=theta, beta=beta, Njs=list(Njs*runif(length(Njs),0.9,1.1))))
#
#estimate.rds(degree.sampled.vec, Sij = make.Sij(degree.sampled.vec), control = generate.rds.control(maxit=10, method="BFGS"),
#		initial.values = list(theta=theta, beta=beta, Njs=list(Njs*2)))







## Grid search given theta and beta parameters:
grid.two.parameters<- function(sampled.degree.vector, Sij, grid.values, Nj.infaltions, control){
	grid.values<- makeNjGrid(sampled.degree.vector, grid.values, Nj.infaltions)
	
	result<- grid.all.parameters(sampled.degree.vector, Sij, grid.values, control)
	return(result)	
}
## Testing:
#Njs<- c(20,20,20,20); names(Njs)<- c("1","50","100","1000"); Njs<- as.table(Njs)
#theta<- 1.2
#beta<- 1.2e-8
#tail(degree.sampled.vec<- generate.sample(theta, Njs, beta, sample.length=1e4))
#plot(degree.sampled.vec, type='h', main='Sample')
#
#beta.deflations<- c(1)
#Nj.infaltions<- c(1,2,10)
#thetas<- theta * c(0.9,1,1.1)
#betas<- beta *c(0.9,1,1.1)
#theta.beta.grid<- makeBetaGrid(degree.sampled.vec, thetas, betas, beta.deflations)
#
#grid.two.parameters(degree.sampled.vec, Sij = make.Sij(degree.sampled.vec),	Nj.infaltions = Nj.infaltions,			
#		grid.values = theta.beta.grid,  
#		control = generate.rds.control( maxit = 20, method="Nelder-Mead"))













## Grid search given theta alone:
grid.one.parameters <- function(sampled.degree.vector, Sij, thetas, beta.inflations=c(1), Nj.inflations, 
		control){
	betas<- c()
	for(theta in thetas) betas<- c(betas, exp(initialLogBetas(sampled.degree.vector, theta)))
	
	theta.beta.values<- makeBetaGrid(sampled.degree.vector, thetas, betas, beta.inflations)
	
	result<- grid.two.parameters(sampled.degree.vector, Sij, theta.beta.values, Nj.infaltions = Nj.inflations, control)
	return(result)
}
## Testing:
#Njs<- c(20,20,20,20); names(Njs)<- c("1","50","100","1000"); Njs<- as.table(Njs)
#theta<- 1.2
#beta<- 1.2e-8
#tail(degree.sampled.vec<- generate.sample(theta, Njs, beta, sample.length=1e4))
#plot(degree.sampled.vec, type='h', main='Sample')
#
#Nj.infaltions<- c(1,2,10)
#thetas<- theta * c(0.9,1,1.1)
#
#grid.one.parameters(degree.sampled.vec, Sij = make.Sij(degree.sampled.vec), 
#		thetas = thetas, Nj.inflations = Nj.infaltions,
#		beta.inflations= c(1, 100, 1e3, 1e5),
#		control = generate.rds.control( maxit = 10, method="Nelder-Mead"))







# Optimize using grid search:
estimate.rds.grid<- function(sampled.degree.vector, Sij, grid.values,		
		control=generate.rds.control(maxit=10, method="Nelder-Mead")){
	
	message('Caution: This might take several hours to run. Consider estimate.rds as a faster (but less robust) alternative.')
	
	Nj.inflations<- control$Nj.inflations
	beta.inflations<- control$beta.inflations
	
	## Initiate using no values: 
	if (missing(grid.values)) {
		grid.values<- seq(-1,1, by=0.1)		
		result<- grid.one.parameters(sampled.degree.vector, Sij = Sij, 
				thetas = grid.values, Nj.inflations = Nj.inflations,
				beta.inflations= beta.inflations,
				control = control)
	}
	
	
	## Grid of theta, beta and Njs:
	else if(length(grid.values[[1]])==3L) {
		result<- grid.all.parameters(sampled.degree.vector, Sij = Sij, grid.values = grid.values, control = control)
	} 	
	

	## Grid given theta and beta
	else if(length(grid.values[[1]])==2L) {
		result<- grid.two.parameters(sampled.degree.vector, Sij = Sij,	
				Nj.infaltions = Nj.inflations,			
				grid.values = grid.values,  
				control = control)
	}
	
				
	
	## Grid given theta alone:
	else if(length(grid.values[[1]])==1L) {
		result<- grid.one.parameters(sampled.degree.vector, Sij = Sij, 
				thetas = grid.values, Nj.inflations = Nj.inflations,
				beta.inflations= beta.inflations,
				control = control)
	}
	
	else stop("Invalid initialization values. See examples.")
	
	
	return(result)	
}
##### Testing with simulated data #### 
#require(chords)
#Njs<- c(100,100,100,100); names(Njs)<- c("1","50","100","1000"); Njs<- as.table(Njs)
#theta<- 1
#beta<- 2.5e-8
#tail(degree.sampled.vec<- generate.sample(theta, Njs, beta, sample.length=1e4))
#x11(); plot(degree.sampled.vec, type='h', main='Sample')
#
## Veryfying likelihood computation:
#estimate.rds(degree.sampled.vec, Sij = make.Sij(degree.sampled.vec), 				
#		initial.values = list(theta=c(theta), beta=c(beta), Njs=list(Njs)),  
#		control = generate.rds.control( maxit = 1, method="BFGS"))
#
##### Grid over (theta,beta,Njs) from true parameter values:
#thetas<- theta * c(0.9,1,1.1)
#betas<- beta *c(0.9,1,1.1)
#theta.beta.grid<- makeBetaGrid(degree.sampled.vec, thetas, betas, 1)
#grid.values<- makeNjGrid(degree.sampled.vec, theta.beta.grid, 1)		
#
#estimate.rds.grid(degree.sampled.vec, Sij = make.Sij(degree.sampled.vec),				
#		grid.values = grid.values, control = generate.rds.control(maxit=200, method="BFGS"))
#
###### Grid over (theta,beta) from true parameter values:
#beta.deflations<- c(1)
#thetas<- theta * c(0.9,1,1.1)
#betas<- beta *c(0.9,1,1.1)
#grid.values<- makeBetaGrid(degree.sampled.vec, thetas, betas, beta.deflations)
#
#estimate.rds.grid(degree.sampled.vec, Sij = make.Sij(degree.sampled.vec),				
#		grid.values = grid.values, control = generate.rds.control(maxit=20, method="Nelder-Mead"))
#
#
#### Grid over (theta) from true parameter values:
#grid.values<- theta * c(0.9,1,1.1)
#
#estimate.rds.grid(degree.sampled.vec, Sij = make.Sij(degree.sampled.vec),				
#		grid.values = grid.values,  
#		control = generate.rds.control( maxit = 10, method="Nelder-Mead"))
#
#### Grid with no starting parameters:
#estimate.rds.grid(degree.sampled.vec, Sij = make.Sij(degree.sampled.vec), control = generate.rds.control( maxit = 2, method="Nelder-Mead"))




#### Testing with true data ####
#require(chords)
#data(simulation, package='chords')
#temp.data<- unlist(data3[1,7000:7500])
## Use Simon's (theta,beta) grid with no Nj inflation:
#estimate.rds.grid(temp.data, Sij = ?!?!!)
















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
