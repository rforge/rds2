rm(list=ls())

#### Generate a network ####
network.size<- 1000
network.density<- 0.02
neighbour.count<- lapply(seq(1, network.size), function(x) rbinom(1, network.size, network.density))
# incidence.list<- lapply(neighbour.count, function(x) sample(x=seq(1,network.size), size=x))

#### Generate a sample under the RDS assumption ####
theta<- 1
Njs<- table(unlist(neighbour.count))
js<- as.numeric(names(Njs))
N<- network.size
(maximal.beta<- 1/ (N * max(Njs) * max(js) ))
(beta<- maximal.beta*200) # Note: larger beta favour sampling non 0 degrees.
## TODO: Find the normalizing constant.

sample.length<- 1000
Uik<- rep(0, length.out=max(js))
names(Uik)<- seq(along.with=Uik)
snowball<- 1
degree.sampled.vec<- rep(NA, sample.length)
for(i in seq(1, sample.length)){
  pi.function<- function(k) beta * k^theta * snowball * (Njs[paste(k)] - Uik[k] )     
  probs.function<- function(k)  pi.function(k) * prod(sapply(js[js!=k], function(x) 1-pi.function(x)  ))
  degree.probabilities<- sapply(js, probs.function)
  (no.sample.probability<- prod(sapply(js, function(x) 1-pi.function(x) )))
  sample.someone<- as.logical(rbinom(1, 1, prob=1-no.sample.probability ))
  
  # if someone is sampled:
  if(sample.someone){
    # sample degree
    degree.sampled<- sample(x=js, prob=degree.probabilities, size=1)      
    # increase Uik
    Uik[degree.sampled]<- Uik[degree.sampled]+1
    snowball<- sum(Uik)
    degree.sampled.vec[i]<- degree.sampled    
  }  
  else{
    degree.sampled.vec[i]<- 0
  }
}
tail(degree.sampled.vec)
head(degree.sampled.vec)


#### Estimation ####
# require(rds2)
# (estimated.simulation<- estimate.rds3(data=sampled.degrees, Sij=make.Sij(data=sampled.degrees), initial.thetas=c(0.5,1,1.5), const=0.5, theta.minimum=-2, theta.range=4))
# comparison(estimated.simulation)
