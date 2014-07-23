source('temp/jass_example.R')
library(dplyr)
rds.sample2 <- tbl_df(rds.sample)

#-------------- Curitiba Data -----------------#

brazil.csv <- read.csv(file='~/Dropbox/RDS/Round2/Brazil/multiplier.csv', stringsAsFactors=FALSE)
brazil.csv2 <- tbl_df(brazil.csv)
brazil.csv3 <- select(brazil.csv2, 
                      uid, netsize.1.bss, InCoupon, OutCoupon1, OutCoupon2, OutCoupon3, interview.date)
names(brazil.csv3)
names(brazil.csv3) <- names(rds.sample2)

library(lubridate)
library(magrittr)
brazil.csv3$interviewDt2 <- brazil.csv3$interviewDt
brazil.csv3$interviewDt <- as.character(brazil.csv3$interviewDt) %>%  dmy %>% unclass

rds.object2<- initializeRdsObject(brazil.csv3)
str(rds.object2)
rds.object2$estimates <- estimate.b.k(rds.object = rds.object2 )
str(rds.object2)
plot(rds.object2$estimates$Nk.estimates, type='h')
plot(rds.object2$estimates$log.bk.estimates, type='h')
getTheta(rds.object2)



