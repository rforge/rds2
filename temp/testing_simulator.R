source('temp/Uganda_example.R')

#---------- Simulate RDS sample --------#

rds.simulated.object <- makeRdsSample(
  N.k =nk.estimates$Nk.estimates , 
  b.k = exp(nk.estimates$log.bk.estimates),
  sample.length = 900L)
str(rds.simulated.object)
plot(rds.simulated.object$degree.in, type='h')
plot(rds.object$rds.sample$NS1, type='h')
plot(rds.simulated.object$rds.sample$interviewDt, cex=0.2)
plot(rds.object$rds.sample$interviewDt, cex=0.2)
nk.estimates.2 <- estimate.b.k(rds.object = rds.simulated.object )
compareNkEstimate(nk.estimates.2, nk.estimates)

sum(nk.estimates$Nk.estimates)
sum(nk.estimates.2$Nk.estimates)


## Repeat and compute average MSE:
replicationss <- 100L
MSEs <- replicate(replicationss,{
  rds.simulated.object <- makeRdsSample(
    N.k =nk.estimates$Nk.estimates , 
    b.k = exp(nk.estimates$log.bk.estimates),
    sample.length = 900L)
  nk.estimates.2 <- estimate.b.k(rds.object = rds.simulated.object )
  MSE <- compareNkEstimate(nk.estimates.2, nk.estimates)
})
dput(MSEs)

