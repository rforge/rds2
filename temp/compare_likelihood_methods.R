rm(list=ls())



#---- Compute on hand made data:---------
N <- 2000# polpulation size
sl <- 1000 # sample size
obs.deg <- 5

rds.sample <- data.frame(MyUniID = seq(1, by = 1, len = sl) ,
                         NS1 =  seq(obs.deg, by = 0, len = sl),
                         refCoupNum  = seq(0, by = 1, len = sl), 
                         coup1 = seq(1, by = 1, len = sl),
                         coup2 = seq(sl+1, by = 1, len = sl), 
                         coup3 = seq(2*sl+1, by = 1, len = sl),
                         interviewDt = seq(0, by = 0, len = sl), 
                         stringsAsFactors = TRUE)
rds.sample$NS1[1] <- 3
for (r in 2:sl) rds.sample$interviewDt[r] <- rds.sample$interviewDt[r-1] + 1/(obs.deg*(r-1)*(N-r+1))

log.bks <- log(1:obs.deg)
Nk.s <- rep(N,obs.deg)
I.t <- 1:sl
n.k.counts <- rep(NA, obs.deg)
.tab <- table(rds.sample$NS1[-1])
n.k.counts[as.numeric(names(.tab))]<- .tab
degree.in <- c(3,rep(obs.deg,sl-1))
degree.out <- rep(0,sl)
arrive.intervals <- diff(rds.sample$interviewDt)
arrival.degree <- degree.in

chords:::likelihood(log.bk = log.bks, 
                    Nk.estimates = Nk.s, 
                    I.t = I.t, 
                    n.k.counts = n.k.counts, 
                    degree.in = degree.in, 
                    degree.out = degree.out, 
                    arrival.intervals = arrive.intervals, 
                    arrival.degree = arrival.degree)


chords:::likelihood.theta.2(
  log.bk = log.bks, 
  Nk.estimates = Nk.s, 
  I.t = I.t, 
  n.k.counts = n.k.counts, 
  degree.in = degree.in, 
  degree.out = degree.out, 
  arrival.intervals = arrive.intervals, 
  arrival.degree = arrival.degree,
  const = 1)





# Compute Jonathan's likelihood at some value:

true.Nks <- rep(0,100); true.Nks[c(2,100)] <- 1000
theta <- 1e-1
true.log.bks <- rep(-Inf, 100)
true.log.bks[c(2,100)] <- theta*log(c(2,100))
sample.length <- 1000L

rds.simulated.object <- makeRdsSample(
  N.k =true.Nks, 
  b.k = exp(true.log.bks),
  sample.length = sample.length)
rds.simulated.object$estimates <- estimate.b.k(rds.simulated.object)

chords:::likelihood(log.bk = rds.simulated.object$estimates$log.bk.estimates, 
                    Nk.estimates = rds.simulated.object$estimates$Nk.estimates, 
                    I.t = rds.simulated.object$I.t, 
                    n.k.counts = rds.simulated.object$estimates$n.k.counts, 
                    degree.in = rds.simulated.object$degree.in, 
                    degree.out = rds.simulated.object$degree.out, 
                    arrival.intervals = rds.simulated.object$estimates$arrival.intervals, 
                    arrival.degree = rds.simulated.object$estimates$arrival.degree)


chords:::likelihood.theta.2(
  log.bk = rds.simulated.object$estimates$log.bk.estimates, 
  Nk.estimates = rds.simulated.object$estimates$Nk.estimates, 
  I.t = rds.simulated.object$I.t, 
  n.k.counts = rds.simulated.object$estimates$n.k.counts, 
  degree.in = rds.simulated.object$degree.in, 
  degree.out = rds.simulated.object$degree.out, 
  arrival.intervals = rds.simulated.object$estimates$arrival.intervals, 
  arrival.degree = rds.simulated.object$estimates$arrival.degree,
  const = 1)


## Try maximum likelihood 
chords:::estimate.b.theta(rds.simulated.object)
