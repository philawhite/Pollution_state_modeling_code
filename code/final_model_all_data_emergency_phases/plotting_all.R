library(MASS)
library(MCMCpack)
library(geoR)
library(dplyr)
library(parallel)
library(xtable)
#library(doMC)
library(foreach)
library(lubridate)

try(setwd("C:/Users/phila/Dropbox/Research/Mexico/code/final3"),silent= TRUE)
source("functions_AR_Final_all.R")
#load("preds.RData")
load("~/Documents/probs_final_all.RData")
#load("~/Documents/preds.RData")

post_sum = function(x,alp = 0.05){
  m = mean(x)
  q = quantile(x,c(alp/2,1-alp/2))
  return(c(m,q))
}


table_funcs = function(phase ,dat_months,dat_dates,dat_hours ,nr ,reps ){
  
  ### nr, nt, reps
  
  nt = dim(phase)[2]
  day_phase = apply(phase,c(1,3),function(x){
    tapply(x,dat_dates,max)
  })
  
  day_mons = month(unique(dat_dates))
  
  temp =  apply(phase,c(2,3),max)
  
  day_temp = apply(temp,2,function(x){
    tapply(x,dat_dates,max)
  })
  
  reg_hour = apply(phase,c(1,3),sum)   ### rows are regions , cols are reps
  reg_day = apply(day_phase,c(2,3),sum)   ### rows are regions , cols are reps
  
  tot_hour = apply(temp,2,sum)
  tot_day = apply(day_temp,2,sum)
  
  reg_hour_mon = apply(phase,c(1,3),function(x){
    tapply(x,dat_months,sum)
  })
  
  reg_day_mon = apply(day_phase,c(2,3),function(x){
    tapply(x,day_mons,sum)
  })
  
  tot_hour_mon = apply(temp,2,function(x){
    tapply(x,dat_months,sum)
  })
  
  tot_day_mon = apply(day_temp,2,function(x){
    tapply(x,day_mons,sum)
  })
  
  reg_hour_hour = apply(phase,c(1,3),function(x){
    tapply(x,dat_hours,sum)
  })
  
  tot_hour_hour = apply(temp,2,function(x){
    tapply(x,dat_hours,sum)
  })
  
  return(list(reg_day = reg_day,reg_hour = reg_hour,tot_day = tot_day,tot_hour = tot_hour, 
              reg_hour_mon = reg_hour_mon,reg_day_mon = reg_day_mon,tot_hour_mon = tot_hour_mon,
              tot_day_mon = tot_day_mon,reg_hour_hour = reg_hour_hour, tot_hour_hour = tot_hour_hour))
}



####################################################
#############   Get summaries  ############
####################################################

dat_dates = as.Date(dat$Date,format = "%m/%d/%Y")[1:(nt*ns) %% 24 ==1]
dat_hours = rep(1:24,365)
dat_months = month(dat_dates)

out_pm10 = table_funcs(phase = exceed_PM10[,-(1:(24*31)),],dat_months = dat_months[-(1:(24*31))],
                       dat_dates= dat_dates[-(1:(24*31))], dat_hours = dat_hours[-(1:(24*31))],
                       nr = nr,reps = reps)  

out_O3 = table_funcs(phase = exceed_O3[,-(1:(24*31)),],dat_months = dat_months[-(1:(24*31))],
                     dat_dates= dat_dates[-(1:(24*31))], dat_hours = dat_hours[-(1:(24*31))],
                     nr = nr,reps = reps)  

####################################################
########################## Day and Hour Summaries
####################################################

tab_cols = c("CE","NE","NW","SE","SW","Any")

PM10_hour = apply(rbind(out_pm10$reg_hour,out_pm10$tot_hour),1,post_sum)
colnames(PM10_hour) = tab_cols

PM10_day = apply(rbind(out_pm10$reg_day,out_pm10$tot_day),1,post_sum)
colnames(PM10_day) = tab_cols



O3_hour = apply(rbind(out_O3$reg_hour,out_O3$tot_hour),1,post_sum)
colnames(O3_hour) = tab_cols

O3_day = apply(rbind(out_O3$reg_day,out_O3$tot_day),1,post_sum)
colnames(O3_day) = tab_cols

xtable(cbind(O3_hour,O3_day),digits = 0)

xtable(cbind(O3_hour/8016,O3_day/334),digits = 4)

xtable(cbind(PM10_hour,PM10_day),digits = 0)

xtable(cbind(PM10_hour/8016,PM10_day/334),digits = 4)

####################################################
########################## Monthly Summaries
####################################################

mon_labs = unique(format(dat_dates, format="%b"))


mon_o3_hour_reg = round(apply(out_O3$reg_hour_mon,c(1,2),post_sum))
mon_o3_hour_tot = round(apply(out_O3$tot_hour_mon,1,post_sum))

day_mon = c(28,31,30,31,30,31,31,30,31,30,31)
hour_mon = 24*day_mon

#mon_o3_day_reg =round(apply(out_O3$reg_day_mon,c(1,2),post_sum))
#mon_o3_day_tot =round(apply(out_O3$tot_day_mon,1,post_sum))


pdf("../../Writing/ozone_bymonth.pdf")
ylims = range(c(apply(mon_o3_hour_reg,3,function(x){range(x[1,])}) ,mon_o3_hour_tot[1,]))
par(mar = c(5,5,3,1))
plot(2:12,mon_o3_hour_reg[,,1][1,]/hour_mon,type="o",col=2,ylim = c(0,1/2),ylab = "Proportion of Hours",xaxt = "n",cex.axis = 1.4,
     main = "Ozone Exceedance by Month",xlab="",cex.lab = 1.8,lwd = 2,cex.main = 1.8,lty = 1)
axis(1,at=2:12, labels=mon_labs[-1],las = 2,col.axis="black",cex.axis=1.8, tck=-.03)

for(i in 2:5) lines(2:12,mon_o3_hour_reg[,,i][1,]/hour_mon,type="o",col=(i+1),lwd = 2,lty = i)
lines(2:12,mon_o3_hour_tot[1,]/hour_mon,col = 1,type="o",lwd = 2,lty = 6)

#legend("topright",c("CE","NE","NW","SE","SW","Any"),col = c(2:6,1),lwd=2,ncol=2,cex=1.4)
dev.off()


mon_pm_hour_reg = round(apply(out_pm10$reg_hour_mon,c(1,2),post_sum))
mon_pm_hour_tot = round(apply(out_pm10$tot_hour_mon,1,post_sum))

#mon_pm_day_reg =round(apply(out_pm10$reg_day_mon,c(1,2),post_sum))
#mon_pm_day_tot =round(apply(out_pm10$tot_day_mon,1,post_sum))


pdf("../../Writing/pm10_bymonth.pdf")
ylims = range(c(apply(mon_pm_hour_reg,3,function(x){range(x[1,])}) ,mon_pm_hour_tot[1,]))
par(mar = c(5,5,3,1))

plot(2:12,mon_pm_hour_reg[,,1][1,]/hour_mon,type="o",col=2,ylim = c(0,1),ylab = "Proportion of Hours",xaxt = "n",cex.axis = 1.4,
     main = expression(paste(PM[10]," Exceedance by Month")),xlab="",cex.lab = 1.8,lwd = 2,cex.main = 1.8,lty = 1)
axis(1,at=2:12, labels=mon_labs[-1],las = 2,col.axis="black",cex.axis=1.8, tck=-.03)

for(i in 2:5) lines(2:12,mon_pm_hour_reg[,,i][1,]/hour_mon,type="o",col=(i+1),lwd = 2,lty = i)
lines(2:12,mon_pm_hour_tot[1,]/hour_mon,col = 1,type="o",lwd = 2,lty = 6)

#legend("top",c("CE","NE","NW","SE","SW","Any"),col = c(2:6,1),lwd=2,ncol=3,cex=1.4)
dev.off()

#################################################################
####################################### Hourly Summaries
#################################################################

hour_labs = paste(c(1:23,0),":00",sep="")


hour_o3_hour_reg = round(apply(out_O3$reg_hour_hour,c(1,2),post_sum))
hour_o3_hour_tot = round(apply(out_O3$tot_hour_hour,1,post_sum))

#hour_o3_day_reg =round(apply(out_O3$reg_day_hour,c(1,2),post_sum))
#hour_o3_day_tot =round(apply(out_O3$tot_day_hour,1,post_sum))


pdf("../../Writing/ozone_byhour.pdf")
ylims = range(c(apply(hour_o3_hour_reg,3,function(x){range(x[1,])}) ,hour_o3_hour_tot[1,]))
par(mar = c(5,5,3,1))

plot(1:24,hour_o3_hour_reg[,,1][1,]/334,type="o",col=2,ylim = c(0,.7),ylab = "Proportion of Hours",xaxt = "n",cex.axis = 1.4,
     main = "Ozone Exceedance by Hour",xlab="",cex.lab = 1.8,lwd = 2,cex.main = 1.8,lty = 1)
axis(1,at=1:24, labels=hour_labs,las = 2,col.axis="black",cex.axis=1.5, tck=-.03)

for(i in 2:5) lines(hour_o3_hour_reg[,,i][1,]/334,type="o",col=(i+1),lwd = 2,lty = i)
lines(hour_o3_hour_tot[1,]/334,col = 1,type="o",lwd = 2,lty = 6)

#legend("topleft",c("CE","NE","NW","SE","SW","Any"),col = c(2:6,1),lwd=2,ncol=2,cex=1.4)
dev.off()

hour_pm_hour_reg = round(apply(out_pm10$reg_hour_hour,c(1,2),post_sum))
hour_pm_hour_tot = round(apply(out_pm10$tot_hour_hour,1,post_sum))

#hour_pm_day_reg =round(apply(out_pm10$reg_day_mon,c(1,2),post_sum))
#hour_pm_day_tot =round(apply(out_pm10$tot_day_mon,1,post_sum))
# pdf("../../Writing/pm10_byhour.pdf")
# ylims = range(c(apply(hour_pm_hour_reg,3,function(x){range(x[1,])}) ,hour_pm_hour_tot[1,]))
# 
# plot(1:24,hour_pm_hour_reg[,,1][1,],type="o",col=2,ylim = c(ylims[1],ylims[2]+30),
#      ylab = "Number of Hours",xaxt = "n",main = expression(paste(PM[10]," Exceedance by Hour")),
#      xlab="",cex.lab = 1.8,lwd = 2,cex.main = 1.8,lty = 1)
# axis(1,at=1:24, labels=hour_labs,las = 2,col.axis="black",cex.axis=1.8, tck=-.03)
# 
# for(i in 2:5) lines(hour_pm_hour_reg[,,i][1,],type="o",col=(i+1),lwd = 2,lty = i)
# lines(hour_pm_hour_tot[1,],col = 1,type="o",lwd = 2,lty = 6)
# 
# legend("top",c("CE","NE","NW","SE","SW","Any"),col = c(2:6,1),lwd=2,ncol=6,cex=.8)
# dev.off()

rm(list = ls())
