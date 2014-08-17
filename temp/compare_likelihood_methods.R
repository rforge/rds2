rm(list=ls())

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
                    arrival.degree = rds.simulated.object$estimates$arrival.degree,
                    const = 20)


## Simon's likelihood:
# dkt = degrees of each individual over time
# tt = timelist
# scale parameters
# logpar = are parameters log transformed? exclude theta
# constrain = if true, n1 and n2 represent the 'excess' over the sampled numbers
likelihood.theta.2 <- function (log.bk, Nk.estimates, I.t, 
                                n.k.counts, degree.in, degree.out, 
                                arrival.intervals, arrival.degree, const=1) {
  Sk <- Nk.estimates[Nk.estimates>0]
  dk <- unique(arrival.degree)
  ndeg <- length(dk)
  N <- sum(Nk.estimates)
  res <- 0
  nstep <- length(arrival.degree)
  # start by subtracting seed
  It <- rep(0,nstep)
  Ut <- rep(0,nstep)
  # choose first obs as seed
  seeddeg <- arrival.degree[1]
  Sk[match(seeddeg,dk)] <- Sk[match(seeddeg,dk)]-1
  It[1] <- 1
  Ut[1] <- 1
  betas <- exp(log.bk[!is.na(log.bk)])
  
  # run through timesteps and calculate loglik
  for(i in 2:nstep){
    dkrates <- betas * It[i-1] * Sk
    # total rates
    sumrates <- sum(dkrates)
    deltat <- arrival.intervals[i-1]
    res <- res + dexp(deltat, sumrates, log=TRUE)
    thisdeg <- arrival.degree[i]
    thisdegidx <- match(thisdeg,dk)
    .prob <- dkrates[thisdegidx]/sumrates
    res <- res + dbinom(1, 1, prob=.prob, log=TRUE)
    Sk[match(thisdeg,dk)] <- Sk[match(thisdeg,dk)]-1
    It[i] <- It[i-1]+1
    Ut[i] <- Ut[i-1]+1
  }
  const*res
}
## Testing:
likelihood.theta.2(
  log.bk = rds.simulated.object$estimates$log.bk.estimates, 
  Nk.estimates = rds.simulated.object$estimates$Nk.estimates, 
  I.t = rds.simulated.object$I.t, 
  n.k.counts = rds.simulated.object$estimates$n.k.counts, 
  degree.in = rds.simulated.object$degree.in, 
  degree.out = rds.simulated.object$degree.out, 
  arrival.intervals = rds.simulated.object$estimates$arrival.intervals, 
  arrival.degree = rds.simulated.object$estimates$arrival.degree,
  const = 1)




#---- Compute on hand made data:---------
N <- 10# polpulation size
sl <- 4 # sample size

rds.sample <- data.frame(MyUniID = seq(1, by = 1, len = sl) ,
                         NS1 =  seq(2, by = 0, len = sl),
                         refCoupNum  = seq(0, by = 1, len = sl), 
                         coup1 = seq(1, by = 1, len = sl),
                         coup2 = seq(sl+1, by = 1, len = sl), 
                         coup3 = seq(2*sl+1, by = 1, len = sl),
                         interviewDt = seq(0, by = 0, len = sl), 
                         stringsAsFactors = TRUE)
rds.sample$NS1[1] <- 3
for (r in 2:sl)
{
  rds.sample$interviewDt[r] <- rds.sample$interviewDt[r-1] + 1/((r-1)*(N-r+1))
}

log.bks <- log(1:3)
Nk.s <- rep(10,3)
I.t <- 1:4
n.k.counts <- c(NA,3)
degree.in <- c(3,rep(2,3))
degree.out <- rep(0,4)
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


likelihood.theta.2(
  log.bk = log.bks, 
  Nk.estimates = Nk.s, 
  I.t = I.t, 
  n.k.counts = n.k.counts, 
  degree.in = degree.in, 
  degree.out = degree.out, 
  arrival.intervals = arrive.intervals, 
  arrival.degree = arrival.degree,
  const = 1)
