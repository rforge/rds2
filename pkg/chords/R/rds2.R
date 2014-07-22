findRecruiter <- function(active.coupons, coupon){
  #     coupon <- 71114
  #      active.coupons <- list()
  #      active.coupons[[1]] <- c(TRUE, TRUE, TRUE)
  #      names(active.coupons[[1]]) <- c(61004,61005,61006)
  #      active.coupons[[2]] <- c(TRUE, TRUE, TRUE)
  #      names(active.coupons[[2]]) <- c(71114,71115,71116)
  
  result <- NA
  recruiter.ind <- sapply(active.coupons, function(x) coupon %in% names(x))
  if(any(recruiter.ind)){
    # Report coupon number:
    recruiter.num <- which(recruiter.ind)
    coupon.ind <- names(active.coupons[[recruiter.num]]) %in% coupon
    coupon.num <- which(coupon.ind)
    
    result <- list(
      recruiter=recruiter.num,
      coupon=coupon.num
    ) 
  }
  return(result)
}
## Testing:
# coupon <- 71114s
# active.coupons <- list()
# active.coupons[[1]] <- c(TRUE, TRUE, TRUE)
# names(active.coupons[[1]]) <- c(61004,61005,61006)
# active.coupons[[2]] <- c(TRUE, TRUE, TRUE)
# names(active.coupons[[2]]) <- c(71114,71115,71116)
# findRecruiter(active.coupons, coupon)






makeSnowBall <- function(rds.sample){
  coupon.inds <- grepl('coup[0-9]*', names(rds.sample))
  sample.length <- ncol(rds.sample)
  
  I.t <- rep(NA, nrow(rds.sample))
  I.t[1] <- 1
  degree.in <- rep(0, nrow(rds.sample))
  degree.out <- rep(0, nrow(rds.sample))
  active.coupons <- list()
  for(period in 2:length(I.t)){
    ### Sketch:
    ## if recruiter in data: 
    # remove incoming coupn
    # update I.t if coupons depleted
    ## if recruiter not in data:
    # update I.t(?)
    ## if coupons handed:
    # add distributed coupons
    # update I.t 
    
    # period <- 2
    # period <- period+1
    I.t.running <- I.t[period-1] 
    
    reference.coupon <- rds.sample[period,'refCoupNum']
    # reference.coupon <- 71114
    recruiter <- findRecruiter(active.coupons, reference.coupon )
    
    if(length(recruiter)>1){
      # revoke coupon
      active.coupons[[recruiter$recruiter]][recruiter$coupon] <- FALSE
      
      # if coupons depleted:
      if(all(!active.coupons[[recruiter$recruiter]])) {
        I.t.running <- I.t.running-1
        degree.out[period] <- attributes(active.coupons[[recruiter$recruiter]])$degree
      }
    }
    
    
    # Increase snowball if coupons handed:
    new.coupons <- rds.sample[period, coupon.inds]
    if(isTRUE(any(new.coupons))){
      new.guy <- rep(TRUE, sum(!is.na(new.coupons)))
      attributes(new.guy) <-list(
        names=new.coupons, 
        degree= rds.sample[period,'NS1'])
      
      active.coupons <- c(active.coupons, list(new.guy))
      
      I.t.running <- I.t.running+1
      degree.in[period] <- rds.sample[period, 'NS1']
    }
    
    I.t[period] <- I.t.running
    
  }
  
  return(list(I.t=I.t,
              degree.in=degree.in,
              degree.out=degree.out))
}
# ## Testing:
# test.snowball <- makeSnowBall(rds.sample)
# str(test.snowball)
# table(rds.sample$NS1)
# table(test.snowball$degree.in)
# table(test.snowball$degree.out)

rdsObjectConstructor <- function(rds.sample=NULL,
                                 I.t=NULL,
                                 degree.in=NULL,
                                 degree.out=NULL,
                                 original.ordering=NULL,
                                 N.k=NULL){
  result <- list(rds.sample=rds.sample,
                 I.t=I.t,
                 degree.in=degree.in,
                 degree.out=degree.out,
                 original.ordering=original.ordering,
                 N.k=N.k)
  return(result)
}
## Testing
# rdsObjectConstructor()

initializeRDS_Object <- function(rds.sample){
  ord <- order(rds.sample[,'interviewDt'])
  rds.sample <-rds.sample[ord,] 
  I.t <- makeSnowBall(rds.sample)
  result <- rdsObjectConstructor(rds.sample=rds.sample,
                                 I.t = I.t$I.t,
                                 degree.in = I.t$degree.in,
                                 degree.out = I.t$degree.out,
                                 original.ordering = ord)
  
  return(result)
}
## Testing:
# rds.object <- initializeRDS_Object(rds.sample)
# ls.str(rds.object)




estimate.b.k.2 <- function(k, A.k, B.k, n.k, n.k.count, k.ind){
  ## Initialize:
  result <- list(N.k=NA)
  N.k <- NA
  
  target <- function(N.k){
    const1 <- N.k-(B.k/A.k)
    pre.const2 <- (N.k - n.k[which(k.ind)-1])
    
    # Deal with impossible N.k:
    if(any(pre.const2<0)) return(-Inf)
    
    const2 <- sum(1/pre.const2)
    const1*const2 - n.k.count
  }
  #   xs <- 0:30; plot(y=sapply(xs, target),x=xs, type='h');abline(0,0)
  
  roots <- NULL
  try(roots <- uniroot(f = target, interval =n.k.count*c(1,1e2)), silent = TRUE)
  if(length(roots)>1) {
    result$N.k <- roots$root
  } else {
    result$N.k <- n.k.count
  }
  
  result$N.k <- as.integer(ceiling(result$N.k))
  return(result)
}
## Testing:
# save(n.k, k.ind,A.k, B.k, file='temp/testing_setup.RData')
# load(file='temp/testing_setup.RData')
# estimate.b.k.2(k = 3, A.k = A.k, B.k = B.k, n.k =n.k, k.ind = k.ind )




## estimate beta_k from sampled degrees and snowball matrix:
estimate.b.k<- function (rds.object) {  
## FIXME: avoid adding k[0] entry that shifts all estimates by 1.
  
  ### Sketch:
  # Generate estimable parameters vector.
  # Optimized parameter-wise.
  
  
  ### Initialize:
  arrival.times <- as.numeric(rds.object$rds.sample$interviewDt)
  ## TODO: adapt import to date formatting.
  arrival.times <- as.integer((arrival.times-min(arrival.times, na.rm = TRUE)) /1e3)
  arrival.intervals <- diff(arrival.times)
  
  arrival.degree<- rds.object$rds.sample$NS1
  max.observed.degree<- max(arrival.degree)
  degree.counts<- table(arrival.degree)
  max.degree.count<- max(degree.counts)
  
  # Sequences per degree
  I.t <- rds.object$I.t
  degree.in <- rds.object$degree.in
  degree.out <- rds.object$degree.out
  
  
  ## Estimate:
  Nk.estimates<- rep(9999L, max.observed.degree) 
  log.bk.estiamtes<- rep(NA, max.observed.degree) 
  names(Nk.estimates)<- max.observed.degree
  uniques<- as.integer(names(degree.counts))
  Nk.estimates[-uniques]<- 0
  for(k in uniques){
    # k <- uniques[[1]]
    k.ind <- arrival.degree==k
    not.k.ind <- !k.ind
    
    # dealing with sample kickoff
    not.k.ind[1] <- FALSE 
    k.ind[1] <- FALSE
    
    n.k <- cumsum((degree.in==k) - (degree.out==k))
    n.k.count <- degree.counts[paste(k)]
    #  tail(cbind(arrival.degree, I.t, n.k, arrival.times))
    #  head(cbind(arrival.degree[k.ind], I.t[which(k.ind)-1], arrival.times[k.ind]))
    #  head(cbind(arrival.degree[k.ind], I.t[which(k.ind)-1], arrival.times[k.ind]))
    
    A.k <- sum( I.t[which(not.k.ind)-1] *  arrival.times[not.k.ind])
    B.k <- sum(n.k[which(not.k.ind)-1] * I.t[which(not.k.ind)-1] *  arrival.times[not.k.ind])    
    
    .temp <- estimate.b.k.2(k=k, A.k=A.k, B.k=B.k, n.k=n.k, 
                            n.k.count= n.k.count, k.ind=k.ind)    
    Nk.estimates[k]<-.temp$N.k
    log.bk.estiamtes[k] <- log(n.k.count)-log(.temp$N.k *  A.k - B.k)
  }
  
  
  result<- list(
    call=sys.call(),
    Nk.estimates=Nk.estimates, 
    log.bk.estimates=log.bk.estiamtes)
  
  return(result)  						
}

## Estimate theta assuming beta_k=beta * k^theta:
getTheta <- function(nk.estimates){
  log.bks <- nk.estimates$log.bk.estimates
  log.ks <- log(seq_along(log.bks))
  lm.1 <- rlm(log.bks~log.ks)
  coefs <- as.list(coef(lm.1) )
  
  return(list(
    log.beta_0 = coefs$`(Intercept)`,
    theta = coefs$log.ks ))
}