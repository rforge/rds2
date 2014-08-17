## Create snowball matrix
makeNKT <- function(uniques, degree.in, degree.out){
  result <- matrix(NA, 
                   nrow=length(uniques), 
                   ncol=length(degree.out), 
                   dimnames=list(k=NULL, t=NULL))
  for(i in seq_along(uniques)){
    #     i <- 36
    result[i,] <- cumsum(degree.in==uniques[i])
  }
  return(result)
}
## Testing:
# makeNKT(uniques, degree.in, degree.out)

likelihoodTheta <- function(
  beta, theta, Nk.estimates, I.t, 
  n.k.counts, degree.in, degree.out, 
  arrival.intervals, arrival.degree){
  ## Verification:
  if(beta < .Machine$double.eps) stop('beta below machine percision.')
  #   if(theta<0) debug()
  ## Computation
  log.beta <- log(beta)
  uniques <- which(!is.na(n.k.counts))
  n.k.t <- makeNKT(uniques, degree.in, degree.out)
  result <- 0.
  for(i in seq_along(arrival.intervals)){
    if(i==1) next()
    for(j in seq_along(uniques)){ 
      #       i <- 5; j <- 5
      #       i <- 5; j <- 5
      k <- uniques[[j]]
      lambda <-  beta * theta^k * (Nk.estimates[k] - n.k.t[j,i-1]) * I.t[i-1]
      lambda <- max(lambda, .Machine$double.eps) 
      
      A <- ifelse(arrival.degree[i]==k, log(lambda), 0) 
      B <- lambda * arrival.intervals[i-1] 
      result <- result + A - B 
    }
  }
  return(result)
}
## Testing:
# example(estimate.b.k)
# theta_0 <- getTheta(rds.object2)
# beta <- exp(theta_0$log.beta_0)
# theta <- theta_0$theta
# chords:::likelihoodTheta(beta, theta, 
#                     rds.object2$estimates$Nk.estimates, 
#                     rds.object2$I.t, 
#                     rds.object2$estimates$n.k.counts, 
#                     rds.object2$degree.in, 
#                     rds.object2$degree.out, 
#                     rds.object2$estimates$arrival.intervals, 
#                     rds.object2$estimates$arrival.degree)






wrap.likelihood <- function(beta, theta, N.k, rds.object){
  I.t <- rds.object$I.t
  n.k.counts <- rds.object$estimates$n.k.counts
  degree.in <- rds.object$degree.in
  degree.out <- rds.object$degree.out
  arrival.intervals <- rds.object$estimates$arrival.intervals
  arrival.degree <- rds.object$estimates$arrival.degree
  
  likelihoodTheta(beta, theta, N.k, I.t, 
             n.k.counts, degree.in, degree.out, arrival.intervals, arrival.degree)
  
}
##Testing:
# likelihoodTheta <- chords:::likelihoodTheta
# chords:::wrap.likelihood(beta, theta, rds.object$estimates$Nk.estimates, rds.object )


estimate.b.theta <-function(rds.object,...){
  ## Initialize:  
  theta_0 <- getTheta(rds.object)
  beta <- exp(theta_0$log.beta_0)
  theta <- theta_0$theta
  N.k <- rds.object$estimates$Nk.estimates
  N.k.ind <- N.k!=0
  
  beta.f <- function(beta) log(beta)
  beta.inv.f <- function(beta.converted) exp(beta.converted )
  
  theta.f <- function(theta) (theta)
  theta.inv.f <- function(theta.converted) (theta.converted)
  
  N.k.f <- function(N.k) log(N.k)
  N.k.inv.f <- function(N.k.converted) exp(N.k.converted)
  
  target <- function(x){
    result <- Inf
    beta <- beta.inv.f(x[1])
    theta <- theta.inv.f(x[2])
    N.k <- rep(0, length(N.k.ind))
    N.k[N.k.ind] <- N.k.inv.f(x[-c(1,2)])
    try(result <- -wrap.likelihood(beta, theta, N.k, rds.object), silent = TRUE)
    return(result)
  }
  
  init <- c(beta.converted=beta.f(theta),
            theta.converted=theta.f(theta), 
            Nks.converted=N.k.f(N.k[N.k.ind]))
# target(init)
  optimal <- optim(par =init ,fn = target, ... )  
  
  new.beta <- beta.inv.f(optimal$par[1])
  new.theta <- theta.inv.f(optimal$par[2])
  new.N.k <- rep(0, length(N.k.ind))
  new.N.k[N.k.ind] <- N.k.inv.f(optimal$par[-c(1,2)])
  #   sum(new.N.k)
  
  result <- list(
    beta=new.beta, 
    theta=new.theta, 
    N.k=new.N.k, 
    optim.result=optimal)
  
  return(result)
} 
## Testing:
# new.estimates <- chords:::estimate.b.theta(rds.object, control=list(maxit=1e3))
# sum(new.estimates$N.k)


