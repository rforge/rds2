rm(list=ls())

#### Generate a network ####
network.size<- 2000
network.density<- 0.02
neighbour.count<- lapply(seq(1, network.size), function(x) rbinom(1, network.size, network.density))
# incidence.list<- lapply(neighbour.count, function(x) sample(x=seq(1,network.size), size=x))

#### Generate a sample under the RDS assumption ####
generate.sample <- function(theta, Njs, beta, sample.length) {
	#sample.length<- 3000
	#theta<- 1
	#Njs<- table(unlist(neighbour.count))
	js<- as.numeric(names(Njs))
	N<- sum(Njs)
	maximal.beta<- 1/ (N * max(Njs) * max(js) )
	#beta<- maximal.beta
	
	error.messge<- paste("beta value too large. Errors might occur. \nFor the given degree frequencies, consider keeping it under ", round(maximal.beta,2), sep="" )
	if(beta > maximal.beta) warning(error.messge) # Note: larger beta favour sampling non 0 degrees.
		
	# Initiate:
	Uik<- rep(0, length.out=max(js))
	names(Uik)<- seq(along.with=Uik)
	snowball<- 1
	degree.sampled.vec<- rep(NA, sample.length)
	# Start sampling:
	for(i in seq(1, sample.length)){	 
		pi.function<- function(k) beta * k^theta * snowball * (Njs[paste(k)] - Uik[k] )
		# Compute the probabilities of sampling degree k=0:
		no.sample.probability<- prod(sapply(js, function(x) 1-pi.function(x) ))
		sample.someone<- as.logical(rbinom(1, 1, prob=1-no.sample.probability ))		
		if(sample.someone){
			# Compute the probabilities of sampling degree k>0:
			probs.function<- function(k)  pi.function(k) * prod(sapply(js[js!=k], function(x) 1-pi.function(x)  ))
			degree.probabilities<- sapply(js, probs.function)
			degree.sampled<- sample(x=js, prob=degree.probabilities, size=1)      
			Uik[degree.sampled]<- Uik[degree.sampled]+1
			snowball<- sum(Uik)
			degree.sampled.vec[i]<- degree.sampled    
		}  
		else{ 
			degree.sampled.vec[i]<- 0
		}
	}
	return(degree.sampled.vec)
}
degree.sampled.vec<-generate.sample(theta = theta, Njs = table(unlist(neighbour.count)), beta = 5e-8, sample.length = 3000)
#save(degree.sampled.vec, file="test.degrees.RData")
#load(file="test.degrees.RData")
plot(degree.sampled.vec, type='h')


#### Estimation ####
 require(rds2)
 (estimated.simulation<- estimate.rds3(data=degree.sampled.vec, Sij=make.Sij(data=degree.sampled.vec), initial.thetas=c(0.5,1,1.5), const=0.5, theta.minimum=-2, theta.range=4))
 # For old version of the package one has to exctract the best result:
best<- estimated.simulation[[which.max(sapply(estimated.simulation, function(x) x$likelihood.optimum))]]
 
 plot(Njs, type='h')	
 points(best$Nj, type='h', col='red')
 
 plot(Njs~best$Nj)
 
 
 
 
