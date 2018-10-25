library(MASS)
library(MCMCpack)
library(geoR)
source("../explore.R")
source("../crps.R")           ### These functions could be replaced by using the package ScoringRules
Rcpp::sourceCpp("Gibbs_updates.cpp")
source("functions_AR.R")

##############################################################################
########################### Model Validation
##############################################################################

reps = 50e3
burn = 5e3
lags = c(1,2,24,7*24)
lags1 <- lags2 <- lags
logY = TRUE
lam = 0
seed = 1

mod = MCMC_mexico(dat = dat,dat_dec = dat_dec,lags1 = lags,lags2 = lags,
                  reps = reps,burn = burn,ns = ns, nt = nt,
                  lam = 0,seed = 1,n_hold = n_hold)

crps1 = numeric(n_hold)
SE1 = numeric(n_hold)
AE1 = numeric(n_hold)
coverage1 = numeric(n_hold)

crps2 = numeric(n_hold)
SE2 = numeric(n_hold)
AE2 = numeric(n_hold)
coverage2 = numeric(n_hold)
ES = numeric(n_hold)


m1 = mean(dat$O3)
m2 = mean(dat$PM10)
s1 = sd(dat$O3)
s2 = sd(dat$PM10)

for(i in 1:n_hold){
  
  crps1[i] = my.crps(mod$hold_dat$O3[i] ,mod$preds1[,i]^2 )
  crps2[i] = my.crps(mod$hold_dat$PM10[i] ,inv_box_tran(mod$preds2[,i],lam1 = 0 ))
  SE1[i] = (mod$hold_dat$O3[i] - mean(mod$preds1[,i]^2))^2
  SE2[i] = (mod$hold_dat$PM10[i] - mean(inv_box_tran(mod$preds2[,i] , lam1 = 0)))^2
  AE1[i] = abs(mod$hold_dat$O3[i] - mean(mod$preds1[,i]^2))
  AE2[i] = abs(mod$hold_dat$PM10[i] - mean(inv_box_tran(mod$preds2[,i],lam1 = 0 )))
  coverage1[i] = in_function(mod$preds1[,i]^2,mod$hold_dat$O3[i])
  coverage2[i] = in_function(inv_box_tran(mod$preds2[,i],lam1 = 0 ),mod$hold_dat$PM10[i] )
  
  ES[i] = my.ES(mod$hold_dat$O3[i],mod$preds1[,i]^2,m1,s1,mod$hold_dat$PM10[i],
                inv_box_tran(mod$preds2[,i],lam1 = 0 ),m2,s2)
}

mean(crps1)
sqrt(mean(SE1))
mean(AE1)
mean(coverage1)

mean(crps2)
sqrt(mean(SE2))
mean(AE2)
mean(coverage2)

mean(ES)


rm(list=setdiff(ls(), c("crps1","crps2","SE1","SE2","AE1","AE2","coverage1","coverage2","ES")))
