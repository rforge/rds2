#### Generate a network ####
rm(list=ls())


network.size<- 1000
network.density<- 0.01
neighbour.count<- lapply(seq(1, network.size), function(x) rbinom(1, network.size, network.density))
#incidence.list<- lapply(neighbour.count, function(x) sample(x=seq(1,network.size), size=x))

#### Generate a sample under the RDS assumption ####
Njs<- table(unlist(neighbour.count[neighbour.count>0]))
theta<- 10
beta<- 1e-19 #2e-8


#### Estimation using package ver 0.65 ####
#install.packages(c("chords_0.65.tar.gz", "~/workspace/rds2/trunk/pkg/chords_0.65.tar.gz"), repos=NULL)
require(chords)
degree.sampled.vec<- generate.sample(theta, Njs, beta, 500000)
plot(degree.sampled.vec, type='h', main='Sample')

(estimated.simulation<- estimate.rds(data=degree.sampled.vec, Sij=make.Sij(data=degree.sampled.vec), initial.thetas=c(0.5,1,1.5), const=20, theta.minimum=-15, theta.range=30))

estimated.simulation$theta
estimated.simulation$beta
cbind(estimated.simulation$Nj, table(degree.sampled.vec))

## Testing the performance ##
pdf('demonstrateEstimates.pdf')
plot(degree.sampled.vec, type='h', main='Sample')
plot(Njs, type='h')	
points(estimated.simulation$Nj, type='h', col='red')
dev.off()

# compare beta*k*theta:
js<- as.numeric(names(Njs))
plot(beta*js^theta~ js, type='h', lwd=2)
points(estimated.simulation$beta * seq(along.with=estimated.simulation$Nj) ^ estimated.simulation$theta, type='h', lwd=1, col='red')











##### Estimation using old package version ####
#require(rds2)
#trace(estimate.rds3, edit='/home/johnros/Applications/Eclipse4.2/eclipse --launcher.openFile') # To edit and trace with Eclipse editor.
#(estimated.simulation<- estimate.rds3(data=degree.sampled.vec, Sij=make.Sij(data=degree.sampled.vec), initial.thetas=c(0.5,1,1.5), const=0.5, theta.minimum=-2, theta.range=4))
##save(estimated.simulation, file="~/Dropbox/Yakir/test.estimation.RData")
##load(file="~/Dropbox/Yakir/test.estimation.RData")
#untrace(estimate.rds3)

 
## Testing the performance:
## Note: The new version of the package dreturns only one component.
#best<- estimated.simulation[[which.max(sapply(estimated.simulation, function(x) x$likelihood.optimum))]]
#pdf('demonstrateEstimates.pdf')
#plot(degree.sampled.vec, type='h', main='Sample')
#plot(Njs, type='h')	
#points(best$Nj, type='h', col='red')
#dev.off()
#
#best$theta
#best$c



 
 
 
 
