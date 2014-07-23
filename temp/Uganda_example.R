#------ Uganda Data (?) -------------#

## Import RDS file with time stamps
col.classes <- c('integer', 'integer', 'integer', 'integer', 'integer', 'integer', 'character')
f.loc <- paste(path.package('chords'), '/extdata/rds_data_Uganda02.csv',sep='')
rds.sample<- read.csv(file = f.loc, colClasses=col.classes)
names(rds.sample)
## Create big rds object:
rds.object<- initializeRdsObject(rds.sample)
str(rds.object)

## Do N.k estimation:
rds.object$estimates <- estimate.b.k(rds.object = rds.object, smooth.degrees.ind = FALSE )
sum(rds.object$estimates$Nk.estimates)
plot(rds.object$estimates$Nk.estimates, type='h')
plot(rds.object$estimates$log.bk.estimates, type='h')
plot(rds.object$estimates$A.ks, type='h')
plot(rds.object$estimates$B.ks, type='h')
plot(rds.object$estimates$n.k.counts, type='h')

## Use b_k=beta*k^theta to smooth population size:
thetaSmoothingNks(rds.object)
