# 
# Author: johnros
###############################################################################














## TODO: B) Pad exported beta_js with NaNs.

estimate.rds.free.betas<- function (sampled.degree.vector, Sij, method="BFGS", initial.values, arc=FALSE, control=generate.rds.control()) {  	
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
		
		N.j<- rep(0, max.observed.degree)
		observed.js.indexes<- 	sapply(Observed.js, function(x) grep(paste("logNjs.",x,"$",sep=""), names(par) )   )
		N.j[Observed.js]<- exp(par[observed.js.indexes]) # fill non trivial Nj estimates.
		beta_js<- exp(par[-observed.js.indexes]) 		
		
		
		# Checking estimates are within allowed range:
		if(isTRUE(   any(round(N.j[Observed.js],10) < round(Observed.Njs,10))  ||   any(beta_js >= 1/(sum(N.j) * N.j[Observed.js]))  ))   {
			return(likelihood.result)
		} 
		
		if(is.numeric(beta_js) &&  all(!is.infinite(beta_js)) ) {
			result<- .C("likelihood_beta_js", 
					sample=as.integer(sampled.degree.vector), 
					Sij=as.integer(Sij),
					S=as.integer(compute.S(Sij)),				  
					beta_js=as.numeric(beta_js), 
					Nj=as.numeric(N.j), 
					observed_degrees=as.integer(rownames(Sij)),
					n=as.integer(length(sampled.degree.vector)), 
					N=as.integer(length(N.j)),
					N_observed=as.integer(nrow(Sij)),				  
					arc=FALSE,
					result=as.double(0))   	  
			likelihood.result<- result$result
		}	  	  
		return(likelihood.result)
	}
	
	
	generate.initial.values<- function(initial.values){
		
		# maximal possible beta:
		initial.log.beta<- log( 1/ (sum(sampled.degree.vector) * Observed.Njs )	)-0.01
		
		wrap.initial.values<- function(...){
			c(
					canonical.betas=initial.log.beta,						
					logNjs=log(Observed.Njs)+0.1
			)
		}
		
		returned.initial.values <- lapply(list(1), wrap.initial.values)		
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
			
			## TODO: B) remove "canonical.beta" names
			result<- list( 
					beta_js= c(exp(x$par[-observed.js.indexes])), 
					initial.values=initial.values,
					Nj=N.j,
					iterations=x$counts,
					likelihood.optimum=x$value, 
					call=the.call)			
		}		
		return(result)
	} 
	
	temp.result<- lapply(likelihood.optim, prepare.result)
	
	length.of.a.proper.output<- 6L
	if(any(sapply(temp.result, length)==length.of.a.proper.output)) {
		clean.temp.result<- temp.result[sapply(temp.result, length)==length.of.a.proper.output]
		final.result<- temp.result[[which.max(sapply(clean.temp.result, function(x) x$likelihood.optimum))]]
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
#(rds.result<- estimate.rds.free.betas(sampled.degree.vector=temp.data , Sij=make.Sij(temp.data), method="BFGS", initial.values=list(), 
#					control=generate.rds.control()))
#str(rds.result)
#plot(rds.result$Nj, type='h', xlab='Degree', ylab=expression(N[j]), main='Estimated Degree Distribution')
#x11()
#plot(rds.result$beta_js ~ sort(unique(temp.data))[-1])
#
#			 
#
### Good values1: network.size = 500, network.density = 0.1, theta<- 3, beta<- 1e-10
#createDegreeCount <- function(network.size, network.density) {
#	neighbour.count<- lapply(seq(1, network.size), function(x) rbinom(1, network.size, network.density))
#	return(table(unlist(neighbour.count[neighbour.count>0])))
#}
#(Njs<- createDegreeCount(network.size = 500, network.density = 0.1))
#theta<- 3
#beta<- 1e-10
#
#
### Manual creation:
##Njs<- c(100,100,100,100); names(Njs)<- c("1","50","100","1000"); Njs<- as.table(Njs)
## Good values2: theta<- 3, beta<- 5e-10
##theta<- 1
##beta<- 5e-8
#
#
#tail(degree.sampled.vec<- generate.sample(theta, Njs, beta, sample.length=1e3))
#plot(degree.sampled.vec, type='h', main='Sample')
#
#
## Testing with initialization from **true** values:
#str(rds.result<- estimate.rds.free.betas(degree.sampled.vec, Sij = make.Sij(degree.sampled.vec), method='BFGS',				
#				initial.values = list(theta=c(theta), beta=c(beta), Njs=list(Njs)),  
#				control = generate.rds.control( maxit = 3000)))
#str(rds.result)
#plot(Njs, type='h', lwd=2.5)	
#points(rds.result$Nj, type='h', col='orange', lwd=2)
#points(rds.result[[1]]$Nj, type='h', col='orange', lwd=2)
#points(rds.result[[2]]$Nj, type='h', col='orange', lwd=2)
#points(rds.result[[3]]$Nj, type='h', col='orange', lwd=2)
#
#js<-sort(unique(degree.sampled.vec))[-1] 
#plot(rds.result$beta_js ~ js)
#
#
#
#
## Initializing with theta only. (true values unknown):
#str(rds.result<- estimate.rds.free.betas(degree.sampled.vec, Sij = make.Sij(degree.sampled.vec), arc=FALSE, 
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
#str(rds.result<- estimate.rds.free.betas(degree.sampled.vec, Sij = make.Sij(degree.sampled.vec), arc=FALSE, 
#				initial.values = list(theta=c(-1, 0.5, 2, 8), beta=c(1e-13, 1e-1, 1e-5)), 
#				method="BFGS", control = generate.rds.control(maxit = 1000)))
#str(rds.result)
#
#
