library(MASS)
library(MCMCpack)
library(geoR)
library(parallel)
#library(doMC)
library(foreach)
library(lubridate)
source("../explore.R")
source("../crps.R")
Rcpp::sourceCpp("Gibbs_updates.cpp")
source("functions_AR_month.R")
file.remove("pred1.txt","pred2.txt")

##############################################################################
########################### Model Validation
##############################################################################

core_use = 4
nr = 5
reps = 10e4
thin = 10
burn = 10e3
lags1 = c(1,2,24,7*24)
lags2 = c(1,2,24,7*24)
mon_pred = 4

##############################################################################
########################### Model Validation
##############################################################################

mod = MCMC_mexico(dat = dat,dat_dec = dat_dec,lags1 = lags1,lags2 = lags2,
          reps = reps, burn = burn,month = mon_pred,core_use = core_use,thin = thin)
reps = reps / thin

rm(list=setdiff(ls(), c("dat","mod","nt","reps","ns","nr","core_use","mon_pred")))
save.image("~/Documents/preds_April.RData")
#save.image("preds.RData")
rm(list=ls())

source("prob_analysis_april.R")
