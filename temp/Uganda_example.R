rm(list=ls())
#------ Uganda Data (?) -------------#

## Import RDS file with time stamps
col.classes <- c('integer', 'integer', 'integer', 'integer', 'integer', 'integer', 'character')
f.loc <- paste(path.package('chords'), '/extdata/rds_data_Uganda02.csv',sep='')
rds.sample<- read.csv(file = f.loc, colClasses=col.classes)
names(rds.sample)
## Create big rds object:
rds.object<- initializeRdsObject(rds.sample)

## Do N.k estimation:
rds.object$estimates <- estimate.b.k(rds.object = rds.object )
sum(rds.object$estimates$Nk.estimates)
plot(table(rds.object$rds.sample$NS1))
plot(rds.object$estimates$Nk.estimates, type='h')

# Are the inter-arrival times exponential?
qs <- length(rds.object$estimates$arrival.intervals) %>% function(x) (1:x)/x
ys <- 1/median(rds.object$estimates$arrival.intervals) %>% qexp(qs, rate = .)
qqplot(x=rds.object$estimates$arrival.intervals, y=ys)


## Use b_k=beta*k^theta to smooth population size:
getTheta(rds.object, robust=TRUE)
sum(thetaSmoothingNks(rds.object, robust=TRUE))






#---- Aggregate degrees to lower variance -----#

bin <- 3
rds.object2<- initializeRdsObject(rds.sample, bin = bin)
rds.object$estimates <- estimate.b.k(rds.object = rds.object2 )
sum(rds.object$estimates$Nk.estimates)

plot(table(rds.object$rds.sample$NS1))
plot(rds.object$estimates$Nk.estimates, type='h')
getTheta(rds.object, bin=bin)




