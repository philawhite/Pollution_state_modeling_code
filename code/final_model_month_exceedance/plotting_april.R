library(MASS)
library(MCMCpack)
library(geoR)
library(dplyr)
library(parallel)
library(xtable)
#library(doMC)
library(foreach)
try(setwd("C:/Users/phila/Dropbox/Research/Mexico/code/final3"),silent= TRUE)
source("functions_AR_month.R")
#load("preds.RData")
load("~/Documents/probs_April.RData")
#load("~/Documents/preds.RData")

###### phase  total over day 
days = as.Date(1:30, origin = "2017-03-31")


pdf("../../Writing/exceed_prob_O3_apr.pdf",width = 10, height = 5)

par(mar = c(5,5,3,1))
plot(days,law_new_prob$reg_O3_day[,1],col=1,type="l",cex.main = 1.5,
     xlab="",xaxt="n",ylim = c(0,1),ylab = "Probabilities",
     cex.lab = 1.8,main = "April Exceedance Probability of Ozone",lty=1,lwd =2)
axis(1,at=as.numeric(days), labels=format(days, format="%m-%d"),
     las = 2,col.axis="blue",
     cex.axis=1.5, tck=-.01)
lines(days,law_new_prob$reg_O3_day[,2],col=2,lty =2,lwd =2)
lines(days,law_new_prob$reg_O3_day[,3],col=3,lty =3,lwd =3)
lines(days,law_new_prob$reg_O3_day[,4],col=4,lty =4,lwd =2)
lines(days,law_new_prob$reg_O3_day[,5],col=5,lty =5,lwd =2)

abline(v = as.numeric(days),lty = 3,col = "lightgray")

dev.off()

pdf("../../Writing/exceed_prob_PM10_apr.pdf",width = 10, height = 5)

par(mar = c(5,5,3,1))
plot(days,law_new_prob$reg_PM10_day[,1],col=1,type="l",cex.main = 1.5,
     xlab="",xaxt="n",ylim = c(0,1),ylab = "Probabilities",
     cex.lab = 1.8,main = "April Exceedance Probability of PM",lty=1,lwd =2)
axis(1,at=as.numeric(days), 
     labels=format(days, format="%m-%d"),las = 2,col.axis="blue",
     cex.axis=1.5, tck=-.01)
lines(days,law_new_prob$reg_PM10_day[,2],col=2,lty =2,lwd =2)
lines(days,law_new_prob$reg_PM10_day[,3],col=3,lty =3,lwd =3)
lines(days,law_new_prob$reg_PM10_day[,4],col=4,lty =4,lwd =2)
lines(days,law_new_prob$reg_PM10_day[,5],col=5,lty =5,lwd =2)

abline(v = as.numeric(days),lty = 3,col = "lightgray")

dev.off()

xx = seq(1/24,30,by = 1/24)

pdf("../../Writing/exceed_prob_O3_apr_hour.pdf",width = 10, height = 5)

par(mar = c(5,5,3,1))
plot(xx,law_new_prob$reg_O3_hour[,1],col=1,type="l",cex.main = 1.5,
     xlab="",xaxt="n",ylim = c(0,1),ylab = "Probabilities",
     cex.lab = 1.8,main = "April Exceedance Probability of Ozone",lty=1,lwd =2)
axis(1,at=seq(1/2,29.5,by = 1)  , 
     labels=format(days, format="%m-%d"),las = 2,col.axis="blue",
     cex.axis=1.5, tck=-.01)
lines(xx,law_new_prob$reg_O3_hour[,2],col=2,lty =2,lwd =2)
lines(xx,law_new_prob$reg_O3_hour[,3],col=3,lty =3,lwd =3)
lines(xx,law_new_prob$reg_O3_hour[,4],col=4,lty =4,lwd =2)
lines(xx,law_new_prob$reg_O3_hour[,5],col=5,lty =5,lwd =2)

abline(v = as.numeric(days),lty = 3,col = "lightgray")

dev.off()

pdf("../../Writing/exceed_prob_PM10_apr_hour.pdf",width = 10, height = 5)

par(mar = c(5,5,3,1))
plot(xx,law_new_prob$reg_PM10_hour[,1],col=1,type="l",cex.main = 1.5,
     xlab="",xaxt="n",ylim = c(0,1),ylab = "Probabilities",
     cex.lab = 1.8,main = "April Exceedance Probability of PM",lty=1,lwd =2)
axis(1,at=seq(1/2,29.5,by = 1)  , 
     labels=format(days, format="%m-%d"),las = 2,col.axis="blue",
     cex.axis=1.5, tck=-.01)
lines(xx,law_new_prob$reg_PM10_hour[,2],col=2,lty =2,lwd =2)
lines(xx,law_new_prob$reg_PM10_hour[,3],col=3,lty =3,lwd =3)
lines(xx,law_new_prob$reg_PM10_hour[,4],col=4,lty =4,lwd =2)
lines(xx,law_new_prob$reg_PM10_hour[,5],col=5,lty =5,lwd =2)

abline(v = as.numeric(days),lty = 3,col = "lightgray")

dev.off()
