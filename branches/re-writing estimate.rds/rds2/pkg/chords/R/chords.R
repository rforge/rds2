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
 



inv.transform.theta<- function(canonical.theta, ...){
	c(canonical.theta, use.names=FALSE)
}


transform.theta<- function(theta, ...){
	c(theta, use.names=FALSE)
}




# Transforms beta back and forth to canonical (real) scale:
generate.beta.scale<- function(scale=10e10) return(scale)

transform.beta <- function(beta, scale=generate.beta.scale()){
	canonical.beta<- log(beta*scale)
	return(c(canonical.beta, use.names=FALSE))
}

inv.transform.beta<- function(canonical.beta, scale=generate.beta.scale()){
	beta<- exp(canonical.beta)/scale
	return(c(beta, use.names=FALSE))
}






# Transforms Nj back and forth to canonical (real) scale:
generate.Nj.scale <- function(scale=1) return(scale)

transform.Nj <- function(Nj, scale=generate.Nj.scale()){
	canonical.Nj<- log(Nj*scale)
	return(canonical.Nj)
}

inv.transform.Nj<- function(canonical.Nj, scale=generate.Nj.scale()){
	Nj<- exp(canonical.Nj)/scale
	return(Nj)
}
## Testing:
#inv.transform.Nj(transform.Nj(100))







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
		fnscale=-1000,
		beta.scale=10e5
		){
	return(list(maxit=maxit, method=method, Nj.inflations=Nj.inflations, beta.inflations=beta.inflations, fnscale=fnscale, beta.scale=beta.scale))
}
## Testing:
#generate.rds.control()




# Generate a list of beta values for each theta and Njs
makeBetas<- function(theta, Nj, Js, beta.inflations){
	result<-list()
	
	beta.bound<- BetaBound(Js, Nj, theta)
	for (inflation in beta.inflations){				
		result<- c(result, list(c(beta=beta.bound*inflation)))			
	}	
	return(result)
}
## Testing:
#makeBetas(theta=1,Nj=c(100,100,100), Js=c(10,50,100), beta.inflations=c(1,1e4,1e8))











makeNjs <- function(sampled.degree.vector, Nj.inflations){
	# Empirical counts
	Observed.Njs <- table(sampled.degree.vector)[-1] # table of degree counts withuot zero
	Observed.js <- as.numeric(names(Observed.Njs))
	
	# Inflate empirical counts (arbitrary)
	result<- list()
	for (inflation in Nj.inflations){
		Njs<-  ceiling(Observed.Njs*inflation)
		names(Njs)<- Observed.js
		Njs.list<- list(Njs)
		result<- c(result, Nj=Njs.list)
	}	
	# Return list of options:
	return(result)
}
### Testing:
#Njs<- c(100,100,100,100); names(Njs)<- c("1","50","100","1000"); Njs<- as.table(Njs)
#theta<- 1
#beta<- 2e-8
#tail(degree.sampled.vec<- generate.sample(theta, Njs, beta, sample.length=1e3))
#makeNjs(degree.sampled.vec, c(1,1.5))












# Generate a list of initialization values:
makeCanonicalInitialization <- function(initial.values, sampled.degree.vector, Nj.inflations, beta.inflations, scale) {
	result<- list()
	
	## Initialize with a thetas only (realistic scenario)	
	if(length(initial.values[[1]])==1L){
		stopifnot(max(sapply(initial.values, length))==min(sapply(initial.values, length)))
		Njs<- makeNjs(sampled.degree.vector, Nj.inflations)		
		# for each theta, compute a grid of betas and Njs					
		for(i in 1:length(initial.values)){
			theta<- initial.values[[i]]$theta										
			for(Nj in Njs){
				betas<- makeBetas(theta=theta, Nj=Nj, Js=as.numeric(names(Nj)), beta.inflations=beta.inflations)
				for(beta in betas){						
					# return in list-of-lists format in canonical scale!
					theta.beta.Nj<- list(
							canonical.theta=transform.theta(theta), 
							canonical.beta=transform.beta(beta,scale), 
							canonical.Nj=transform.Nj(Nj))					
					result<- c(result, list(theta.beta.Nj))						
				}					
			}			
		}			
	}
	
	
	## Initialize with theta and beta (for grid searching)		
	else if(length(initial.values[[1]])==2L){
		stopifnot(max(sapply(initial.values, length))==min(sapply(initial.values, length)))
		Njs<- makeNjs(sampled.degree.vector, Nj.inflations)	
		
		# for each theta and beta, compute a grid of betas and Njs				
		for(i in 1:length(initial.values)){
			theta<- initial.values[[i]]$theta
			beta<- initial.values[[i]]$beta
			for(Nj in Njs){					
				# return in list-of-lists format in canonical scale!
				theta.beta.Nj<- list(
						canonical.theta=transform.theta(theta), 
						canonical.beta=transform.beta(beta, scale), 
						canonical.Nj=transform.Nj(Nj))					
				result<- c(result, list(theta.beta.Nj))						
			}			
		}			
	}
	
	## Intialize with all values (for simulation verification)
	else if(length(initial.values[[1]])==3L){
		stopifnot(max(sapply(initial.values, length))==min(sapply(initial.values, length)))
		# for each theta and beta, compute a grid of betas and Njs				
		for(i in 1:length(initial.values)){
			theta<- initial.values[[i]]$theta
			beta<- initial.values[[i]]$beta
			Nj<- initial.values[[i]]$Njs				
			# return in list-of-lists format in canonical scale!
			theta.beta.Nj<- list(
					canonical.theta=transform.theta(theta), 
					canonical.beta=transform.beta(beta, scale), 
					canonical.Nj=transform.Nj(Nj))					
			result<- c(result, list(theta.beta.Nj))						
		}			
	}
	
	else stop("Bad initialization values")
	
	return(result)
}
### Testing:
#Njs<- c(100,100,100,100); names(Njs)<- c("1","50","100","1000"); Njs<- as.table(Njs)
#theta<- 1
#beta<- 2e-8
#tail(degree.sampled.vec<- generate.sample(theta, Njs, beta, sample.length=1e3))
#makeCanonicalInitialization(list(list(theta=10), list(theta=20)), degree.sampled.vec, c(1,2,3), c(1,2), 1)[[2]]
#
#makeCanonicalInitialization(list(list(theta=10)), degree.sampled.vec, 1, 1, 1e15)[[1]]
#makeCanonicalInitialization(list(list(theta=10, beta=2)), degree.sampled.vec, c(1,2,3), c(1,2), 1)[[1]]
#Njs<-rpois(20,10)
#names(Njs)<- seq(along.with=Njs)
#makeCanonicalInitialization(list(list(theta=10, beta=2, Njs=Njs)), degree.sampled.vec, c(1,2,3), c(1,2), 1)[[1]]
	






 











			 

createDegreeCount <- function(network.size, network.density) {
	neighbour.count<- lapply(seq(1, network.size), function(x) rbinom(1, network.size, network.density))
	return(table(unlist(neighbour.count[neighbour.count>0])))
}
## Testing: 
#(Njs<- createDegreeCount(network.size = 1e5, network.density = 0.01))







# Creates loose bound on maximal (log) beta values given theta and data:
BetaBound <- function(Js, Nj, theta) {
	# Njs assumed to be a table with degree (j) and counts (Njs)
	max.observed.degree<- max(Js)
	maximal.degree.count<- max(Nj)
	pop.size<- Nj %*% Js
	beta<- 1/ (pop.size * maximal.degree.count * max.observed.degree^theta )
	return(beta)
}
## Testing:
#BetaBound(sampled.degree.vector = rpois(200,2), theta = 1)










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
	Uik<- rep(0L, length.out=max(js))
	names(Uik)<- seq(along.with=Uik)
	snowball<- 1L
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
			Uik[degree.sampled]<- Uik[degree.sampled]+1L
			snowball<- sum(Uik)
			degree.sampled.vec[i]<- as.integer(degree.sampled)    
		}  
		else{ 
			degree.sampled.vec[i]<- 0L
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


# Solve for a fixed Nj given theta and beta:
NjSolve <- function(sampled.degree.vector, S, Sij,j, Nj.table, beta, theta, maximal.Nj=1e6, ...){
	const1 <- beta * S[j] * j^theta
	target <- function(Nj){		
		Uij <- Nj-Sij[paste(j),]
		sum( (sampled.degree.vector==j) * 1/Uij - (sampled.degree.vector!=j) * const1 / (1 - const1 * Uij) )
	}
	
	result <- uniroot(target, interval = c(Nj.table[paste(j)]+1, maximal.Nj),...)
	
	return(result)
}
## Testing:
Njs<- c(100,500,500,200)
names(Njs)<- c("1","50","100","1000")
Njs<- as.table(Njs)
theta<- 1.1
beta<- 3e-8
tail(degree.sampled.vec<- generate.sample(theta, Njs, beta, sample.length=1e3))
(.Nj.table <- table(degree.sampled.vec))

matplot(t(.Sij <- make.Sij(degree.sampled.vec)), type="s")
plot(.S <- compute.S(.Sij), type="s")


NjSolve(degree.sampled.vec, .S, .Sij, Nj.table = .Nj.table, j=50, beta = beta, theta = theta, maximal.Nj = 1e10)


