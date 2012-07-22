rm(list=ls())

#### Generate a network ####
network.size<- 1000
network.density<- 0.02
neighbour.count<- lapply(seq(1, network.size), function(x) rbinom(1, network.size, network.density))
incidence.list<- lapply(neighbour.count, function(x) sample(x=seq(1,network.size), size=x))
head(incidence.list)




#### Generate a sample under the RDS assumption ####
sample.length<- 10000
theta<- 1
beta<- 1e-10
Njs<- table(unlist(neighbour.count))
  
sampled<- rep(FALSE, network.size)
sampled[sample(seq(1, network.size),size=1)]<- TRUE
sampled.degrees<- rep(0, sample.length)
for(i in seq(1, sample.length)){
  # compute the probability of each individual being sampled  
  sample.indexes<- seq(1, network.size)[!sampled]
  snowball<- sum(sampled)
  sampled.table<- table(unlist(neighbour.count[sampled]))  
  probs.function<- function(k) {
    # How many subject with degree k have been sampled?
    Uik<- sampled.table[paste(k)]
    Uik<- ifelse(is.na(Uik), 0, Uik)
    # How many subject with degree k are out there?
    Sik<- Njs[paste(k)]-Uik
    return(beta * snowball * Sik * k^theta )
  }  
  probs.vector<- sapply(neighbour.count[!sampled], probs.function)
  
  # should someone be sampled?
  someone.sampled<- as.logical(rbinom(n=1, size=1, prob=1-sum(probs.vector)))  
    # sample with those probabilities
  if(someone.sampled){
    next.sample<- sample(sample.indexes, 1, prob=probs.vector)
    sampled[next.sample]<- TRUE
    sampled.degrees[i]<- neighbour.count[[next.sample]]    
    sample.indexes<- seq(1, network.size)[!sampled]
  }
  else {
    sampled.degrees[i]<- 0
  }  
}

require(rds2)
(estimated.simulation<- estimate.rds3(data=sampled.degrees, Sij=make.Sij(data=sampled.degrees), initial.thetas=c(0.5,1,1.5), const=0.5, theta.minimum=-2, theta.range=4))
comparison(estimated.simulation)



