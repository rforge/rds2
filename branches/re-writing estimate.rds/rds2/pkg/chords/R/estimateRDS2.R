# 
# Author: johnros
###############################################################################


## Estimate using a grid over theta and beta and gradient over Nj:
estimate.rds2<- function (sampled.degree.vector, Sij, initial.values, arc=FALSE, control=generate.rds.control(), all.solutions=FALSE) {  	
	# Initializing:
	the.call<- sys.call()
	method <- control$method
	maxit <- control$maxit
	fnscale <- control$fnscale
	Nj.inflations<- control$Nj.inflations
	beta.inflations<- control$beta.inflations
	beta.scale <- control$beta.scale
	
	max.observed.degree<- max(sampled.degree.vector)
	N.j<- rep(0, max.observed.degree)
	# Look for the degrees in the data with non trivial estimates:
	Observed.Njs<-table(sampled.degree.vector)[-1] # table of degree counts withuot zero
	Observed.js<- as.numeric(names(Observed.Njs))
	maximal.degree.count<- max(Observed.Njs)
	final.result<- NA
	
	
	# Compute the size of the snowball along the sample:
	S<- compute.S(Sij)
	
	
	# Wrap the likelihood function. 
	# Implements constraints on parameters by taking real valued (canonical) parameters and remapping them.
	likelihood.wrap<- function(par){
		# Initialize:
		likelihood.result<- -9999999
		
		beta<- inv.transform.beta(par[['canonical.beta']]) # assumes beta is non negative and given in log scale
		theta<- inv.transform.theta(canonical.theta = par[['canonical.theta']], control)
		
		N.j<- rep(0, max.observed.degree)
		get.j.index<- function(x) grep(paste("canonical.Nj.",x,"$",sep=""), names(par) ) # Get the indexes in par of the Njs.  
		observed.js.indexes<- 	sapply( Observed.js,  get.j.index)
		N.j[Observed.js]<- inv.transform.Nj(par[observed.js.indexes]) # fill non trivial Nj estimates.
		
		
		# Checking estimates are within allowed range:
		bad.Nj<- any( round(N.j[Observed.js],10) < round(Observed.Njs,10) ) # Estimated Nj smaller than observed
		bad.beta<- beta >= BetaBound(Js= Observed.js, Nj=N.j[Observed.js], theta=theta)  # Defined beta values
		if(isTRUE( bad.Nj || bad.beta)) return(likelihood.result) 
		
		# Compute likelihood:
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
	
	
	### Generate initial values:
	# Fills and formats initialization values as required by optim.wrap()
	# Assumes initial.values is a list of lists. Each containing a theta, beta and an Nj vector.	
	
	initial.value.list <- makeCanonicalInitialization(
			initial.values = initial.values, 
			sampled.degree.vector = sampled.degree.vector, 
			Nj.inflations = Nj.inflations, 
			beta.inflations = beta.inflations, 
			scale = beta.scale)
	
	# In case of convergence problems, try different optimization methods. In particular: "Nelder-Mead", "SANN",
	optim.wrap<- function(x) {
		try(optim(par=unlist(x), fn=likelihood.wrap, method=method , control=list(fnscale=fnscale, maxit=maxit)))
	}
	
	likelihood.optim<-lapply(initial.value.list,  optim.wrap )
	
	
	# Convert output from canonical form to original form:
	prepare.result<- function(x){
		result<- simpleError("Optim did not converge")
		
		if(length(x)==5L){
			observed.js.indexes<- unlist(lapply(Observed.js, function(y) grep(paste("canonical.Nj.",y,"$",sep=""), names(x$par))))
			N.j[Observed.js]<- exp(x$par[observed.js.indexes]) # fill non trivial Nj estimates.		
			
			result<- list( 
					beta=inv.transform.beta(x$par[['canonical.beta']]), 
					theta=inv.transform.theta(x$par[['canonical.theta']], control),
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
		stop('Estimation did not converge. Try differet intialization values or optimization method.')
	}
	
	return(final.result)							
}

### Testing: 
## Testing with initialization from **true** values:
# Generate data:
Njs<- c(100,100,100,100)
names(Njs)<- c("1","50","100","1000")
Njs<- as.table(Njs)
theta<- 1.1
beta<- 2.5e-8
tail(degree.sampled.vec<- generate.sample(theta, Njs, beta, sample.length=1e3))
plot(degree.sampled.vec, type='h', main='Sample')

# Estimate:


makeCanonicalInitialization(
		initial.values = lapply(seq(-1,1,0.1), function(x) {list(theta = x)}), 
		degree.sampled.vec, 
		beta.inflations = 1, Nj.inflations = 1, scale=1)

str(rds.result<- estimate.rds2(degree.sampled.vec, Sij = make.Sij(degree.sampled.vec), 
				method='BFGS',				
				initial.values = list(),  
				control = generate.rds.control( maxit = 100)))

# Inspect estimation:
str(rds.result)
plot(Njs, type='h', lwd=2.5)	
points(rds.result$Nj, type='h', col='orange', lwd=2)
points(rds.result[[1]]$Nj, type='h', col='orange', lwd=2)
points(rds.result[[2]]$Nj, type='h', col='orange', lwd=2)
points(rds.result[[3]]$Nj, type='h', col='orange', lwd=2)








