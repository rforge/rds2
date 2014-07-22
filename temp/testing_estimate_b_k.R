#------ Jazz Data (?) -------------#

## Import RDS file with time stamps
col.classes <- c('integer', 'integer', 'integer', 'integer', 'integer', 'integer', 'character')
rds.sample<- read.csv(file = 'inst/extdata/rds_data_Uganda02.csv', colClasses=col.classes)

## Create big rds object:
rds.object<- initializeRDS_Object(rds.sample)

## Do N.k estimation:
nk.estimates <- estimate.b.k(rds.object = rds.object )
str(rds.object)
plot(nk.estimates$Nk.estimates, type='h')
plot(nk.estimates$log.bk.estimates, type='h')

## Testing:
getTheta(nk.estimates)






#---------- Simulate RDS sample --------#

rds.simulated.object <- makeRdsSample(
  N.k =nk.estimates$Nk.estimates , 
  b.k = exp(nk.estimates$log.bk.estimates),
  sample.length = 900L)
plot(rds.simulated.object$degree.in, type='h')
plot(rds.object$rds.sample$NS1, type='h')
plot(rds.simulated.object$rds.sample$interviewDt, cex=0.2)
plot(rds.object$rds.sample$interviewDt, cex=0.2)
nk.estimates.2 <- estimate.b.k(rds.object = rds.simulated.object )
compareNkEstimate(nk.estimates.2, nk.estimates)

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





#-------------- Curitiba Data -----------------#