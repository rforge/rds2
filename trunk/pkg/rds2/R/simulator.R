#### Generate a network ####
rm(list=ls())


network.size<- 2000
network.density<- 0.002
neighbour.count<- lapply(seq(1, network.size), function(x) rbinom(1, network.size, network.density))
#incidence.list<- lapply(neighbour.count, function(x) sample(x=seq(1,network.size), size=x))

#### Generate a sample under the RDS assumption ####
generate.sample <- function(theta, Njs, beta, sample.length) {
	#sample.length<- 3000
	#theta<- 1
	#Njs<- table(unlist(neighbour.count))
	js<- as.numeric(names(Njs))
	N<- sum(Njs)
	maximal.beta<- 1/ (N * max(Njs) * max(js)^theta )
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
		sample.someone<- as.logical(rbinom(1, 1, prob = 1 - no.sample.probability ))		
		if(sample.someone){
			# Compute the probabilities of sampling degree k>0:
			probs.function<- function(k)  pi.function(k) * prod(sapply(js[js!=k], function(x) 1-pi.function(x)  ))
			degree.probabilities<- sapply(js, probs.function)
      #capture.output(sum(degree.probabilities)+ no.sample.probability, file='Yakir.txt', append=TRUE)
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

Njs<- table(unlist(neighbour.count))
theta<- 1
beta<- 9e-9 #2e-8

degree.sampled.vec<-generate.sample(theta = theta, Njs = Njs, beta = beta, sample.length = 10000)
plot(degree.sampled.vec, type='h')
sample.ind<- degree.sampled.vec>0
sample.indexes<- seq(along.with=sample.ind)[sample.ind]
lines(lowess(degree.sampled.vec[sample.ind]~sample.indexes), col='red')
#save(degree.sampled.vec, file="~/Dropbox/Yakir/test.degrees.RData")
#load(file="~/Dropbox/Yakir/test.degrees.RData")



#### Estimation ####
require(rds2)
(estimated.simulation<- estimate.rds3(data=degree.sampled.vec, Sij=make.Sij(data=degree.sampled.vec), initial.thetas=c(0.5,1,1.5), const=0.5, theta.minimum=-2, theta.range=4))
#save(estimated.simulation, file="~/Dropbox/Yakir/test.estimation.RData")
#load(file="~/Dropbox/Yakir/test.estimation.RData")
 
# Testing the performance:
# Note: For old version of the package one has to exctract the best result:
best<- estimated.simulation[[which.max(sapply(estimated.simulation, function(x) x$likelihood.optimum))]]
plot(Njs, type='h')	
points(best$Nj, type='h', col='red')

comparison.1<- cbind(true= beta * as.numeric(names(Njs))^theta, estimated=best$c * best$Nj[as.numeric(names(Njs))]^best$theta )
matplot(comparison.1, type='o')	


 
 
 
 
