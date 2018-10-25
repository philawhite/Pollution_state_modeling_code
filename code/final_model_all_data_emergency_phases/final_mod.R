library(MASS)
library(MCMCpack)
library(geoR)
library(parallel)
#library(doMC)
library(foreach)

source("../explore.R")
source("../crps.R")  ### These functions could be replaced by using the package ScoringRules
Rcpp::sourceCpp("Gibbs_updates.cpp")
source("functions_AR_Final.R")
file.remove("pred1.txt","pred2.txt")

##############################################################################
##############################################################################

core_use = 12
nr = 5
reps = 10e4
thin = 10
burn = 10e3
lags1 = c(1,2,24,7*24)
lags2 = c(1,2,24,7*24)

##############################################################################
########################### Model Validation
##############################################################################

mod = MCMC_mexico(dat = dat,dat_dec = dat_dec,lags1 = lags1,lags2 = lags2,
          reps = reps, burn = burn,ns = ns,nt = nt,core_use = core_use,thin = thin)

reps = reps / thin

rm(list=setdiff(ls(), c("dat","mod","nt","reps","ns","nr","core_use")))
save.image("~/Documents/preds_final.RData")
#save.image("preds.RData")
rm(list=ls())

source("prob_analysis.R")
