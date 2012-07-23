rm(list=ls())

#### Generate a network ####
network.size<- 1000
network.density<- 0.02
neighbour.count<- lapply(seq(1, network.size), function(x) rbinom(1, network.size, network.density))
incidence.list<- lapply(neighbour.count, function(x) sample(x=seq(1,network.size), size=x))
head(incidence.list)




#### Generate a sample under the RDS assumption ####
sample.length<- 2000
theta<- 1
Njs<- table(unlist(neighbour.count))
N<- network.size
js<- as.numeric(names(Njs))

(maximal.beta<- 1/ (N * max(Njs) * max(js) ))
beta<- maximal.beta*20

Uik<- rep(0, length.out=max(js))
names(Uik)<- seq(along.with=Uik)

snowball<- 1
for(i in seq(1, sample.length)){
  # i<- 10
  # compute the probability of degree 0.
  pi.function<- function(k) {
    beta * k^theta * snowball * (Njs[paste(k)]-Uik[k] )    
  }
  
  probs.function<- function(k){
    pi.function(k) * prod(sapply(js[js!=k], pi.function  ))
  }
  
  (no.sample.probability<- 1 - sum(sapply(js, probs.function)))
  
  
  
  # if someone is sampled:
  sample.someon<- as.logical(rbinom(1, 1, prob=no.sample.probability ))
  if(sample.someone){
    Uik[]
    snowball<- snowball+1
  }
  # compute the probability of each degree    
  # update the table of degree counts (Njs)
  snowball<- sum(Uik)  
}
tail(sampled.degrees)
head(sampled.degrees)


#### Estimation ####
require(rds2)
(estimated.simulation<- estimate.rds3(data=sampled.degrees, Sij=make.Sij(data=sampled.degrees), initial.thetas=c(0.5,1,1.5), const=0.5, theta.minimum=-2, theta.range=4))
comparison(estimated.simulation)








# Code from previous versions:
probs.function<- function(k) {
  # How many subject with degree k have been sampled?
  Uik<- sampled.table[paste(k)]
  Uik<- ifelse(is.na(Uik), 0, Uik)
  # How many subject with degree k are out there?
  Sik<- Njs[paste(k)]-Uik
  return( snowball * Sik * k^theta )
}