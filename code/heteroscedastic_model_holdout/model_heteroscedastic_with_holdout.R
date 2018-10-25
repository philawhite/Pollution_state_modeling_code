library(MASS)
library(MCMCpack)
library(geoR)
source("../explore.R")
source("../crps.R")  ### These functions could be replaced by using the package ScoringRules
Rcpp::sourceCpp("Gibbs_updates_het.cpp")
source("functions_AR_daily_heteroscedastic.R")

##############################################################################
########################### Model Validation
##############################################################################

reps = 50e3
burn = 5e3
lags = c(1,2,24,24*7)
lags1 <- lags2 <- lags
logY = TRUE
lam = 0
seed = 1
nh = 24

mod = MCMC_mexico(dat = dat,dat_dec = dat_dec,lags1 = lags,
                  lags2 = lags,reps = reps,burn = burn,ns = ns, nt = nt,nh = nh,logY = FALSE,
                  lam = 0,seed = 1,n_hold = n_hold)

crps1 = numeric(n_hold)
SE1 = numeric(n_hold)
AE1 = numeric(n_hold)
coverage1 = numeric(n_hold)

crps2 = numeric(n_hold)
SE2 = numeric(n_hold)
AE2 = numeric(n_hold)
coverage2 = numeric(n_hold)

for(i in 1:n_hold){
  
  crps1[i] = my.crps(mod$hold_dat$O3[i] ,mod$preds1[,i])
  crps2[i] = my.crps(mod$hold_dat$PM10[i] ,mod$preds2[,i])
  SE1[i] = (mod$hold_dat$O3[i] - mean(mod$preds1[,i]))^2
  SE2[i] = (mod$hold_dat$PM10[i] - mean(mod$preds2[,i]))^2
  AE1[i] = abs(mod$hold_dat$O3[i] - mean(mod$preds1[,i]))
  AE2[i] = abs(mod$hold_dat$PM10[i] - mean(mod$preds2[,i]))
  coverage1[i] = in_function(mod$preds1[,i],mod$hold_dat$O3[i] )
  coverage2[i] = in_function(mod$preds2[,i],mod$hold_dat$PM10[i] )
  
}

mean(crps1)
sqrt(mean(SE1))
mean(AE1)
mean(coverage1)

mean(crps2)
sqrt(mean(SE2))
mean(AE2)
mean(coverage2)


rm(list=setdiff(ls(), c("crps1","crps2","SE1","SE2","AE1","AE2","coverage1","coverage2")))