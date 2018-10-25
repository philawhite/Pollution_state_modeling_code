library(MASS)
library(MCMCpack)
library(geoR)
library(dplyr)
library(parallel)
#library(doMC)
library(foreach)
source("functions_AR_month.R")
#load("preds.RData")
load("~/Documents/preds_April.RData")

# rm(list=ls())
# load("~/Documents/probs_April.RData")
# source("functions_AR_month.R")


phase_all = mclapply(1:reps,function(i){
  Phase_MC(dat = dat,preds = mod$preds[[i]],ns = ns,month = mon_pred,L1_o= 154,L2_o = 204,L1_p = 214,L2_p = 354)
},mc.cores = core_use)

dat_dates = as.Date(dat$Date,format = "%m/%d/%Y")
mon = month(dat_dates)
idx_test = which(mon == mon_pred)  
nt_use = length(unique(dat$time_ind[idx_test]))

mc_phases = array(0,c(nr,nt_use,reps))
exceed_PM10 = array(0,c(nr,nt_use,reps)) 
exceed_O3 = array(0,c(nr,nt_use,reps))

for(i in 1:reps){
  mc_phases[,,i] = t(phase_all[[i]]$phase)
  exceed_PM10[,,i] = t(phase_all[[i]]$exceed_PM10)
  exceed_O3[,,i] = t(phase_all[[i]]$exceed_O3)
}

phase_prob = prob_phase(phase = mc_phases,nr = nr,reps = reps)

law_new_prob = prob_exceed(phase_O3 = exceed_O3,phase_PM10 = exceed_PM10,
                           nr = nr, reps = reps)

save.image("~/Documents/probs_April.RData")
rm(list=ls())

