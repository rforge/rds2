rm(list=ls())
source('temp/Uganda_example.R')
rds.object$estimates <- estimate.b.k(rds.object = rds.object, impute.Nks = FALSE )

#---------- Simulate RDS sample --------#
scaler1 <- 10
rds.object2 <- rds.object
rds.object2$estimates$Nk.estimates <- rds.object$estimates$Nk.estimates* scaler1
rds.simulated.object <- makeRdsSample(
  N.k =rds.object$estimates$Nk.estimates , 
  b.k = exp(rds.object$estimates$log.bk.estimates),
  sample.length = 1e3)
rds.simulated.object$estimates <- estimate.b.k(rds.object = rds.simulated.object )
sum(rds.object2$estimates$Nk.estimates)
sum(rds.simulated.object$estimates$Nk.estimates)
table(rds.simulated.object$estimates$convergence)
chords:::compareNkEstimate(rds.simulated.object$estimates, rds.object$estimates)


sum(rds.object$estimates$Nk.estimates)
sum(rds.simulated.object$estimates$Nk.estimates)



## Repeat and compute average MSE:


replicationss <- 10L
MSEs <- replicate(replicationss,{
  rds.simulated.object <- makeRdsSample(
    N.k =rds.object$estimates$Nk.estimates , 
    b.k = exp(rds.object$estimates$log.bk.estimates),
    sample.length = 1000L)
  rds.simulated.object$estimates <- estimate.b.k(rds.object = rds.simulated.object )
  chords:::compareNkEstimate(rds.simulated.object$estimates, rds.object$estimates)
})
MSEs

