library(rgdal)
library(fields)
#library(forecast)
rm(list=ls())

hd <- function(long1, lat1, long2, lat2) {
  R <- 6371 # Earth mean radius in km - I should change to something else
  d_long <- (long2 - long1)*pi/180
  d_lat <- (lat2 - lat1)*pi/180
  a <- sin(d_lat/2)^2 + cos(lat1*pi/180) * cos(lat2*pi/180) * sin(d_long/2)^2
  c <- 2 * atan2(sqrt(a),sqrt(1-a))
  d <- R * c
  return(d) # Distance in km
}

###### Make the Coordinates

# coords = read.csv("../../data/UTM_Coordinates.csv")
# temp <- data.frame(lon =coords[,2], lat =coords[,3])
# coordinates(temp) = c("lon","lat")
# proj4string(temp) = CRS("+proj=utm +zone=14 ellps=WGS84") #UTM zone 14
# CRS.new = CRS("+init=epsg:4326") # WGS84
# out = spTransform(temp,CRS.new)
# coords = cbind(coords,out@coords)
# rm(CRS.new,out,temp)
# ns = nrow(coords)
# 
# save.image("coord_raw.RData")

###### Load in all the other data
load("coord_raw.RData")
O3_mat = read.csv("../../data/O3.csv")
PM10_mat = read.csv("../../data/PM10.csv")
RH_mat = read.csv("../../data/RH.csv")
TMP_mat = data.matrix(read.csv("../../data/TMP.csv"))
nt = nrow(O3_mat)

O3_mat_dec = read.csv("../../data/O3_dec.csv")
PM10_mat_dec = read.csv("../../data/PM10_dec.csv")

########## Create a single data

dat = matrix(ncol=13,nrow= ns*nt)
colnames(dat) = c("obs_ind","Date","Hour","time_ind","Station","Region","loc_ind","Lon","Lat","PM10","O3","RH","TMP")

obs_ind = 1:(ns*nt)
Date = O3_mat[rep(1:nt,each=ns),"FECHA"]
Hour = O3_mat[rep(1:nt,each=ns),"HORA"]
time_ind = rep(1:nt,each= ns)
Station = rep(coords[,"Stations"],times=nt)
Region = rep(coords[,"Region"],times=nt)
Lat = rep(coords[,"lat"],times=nt)
Lon = rep(coords[,"lon"],times=nt)
loc_ind = rep(1:ns,times= nt)
PM10 <- O3 <- RH <- TMP <- numeric(nt*ns)
count = 1

for(j in 1:nt){
  for(i in 1:ns){
    PM10[count] = PM10_mat[j,2+i]
    O3[count] = O3_mat[j,2+i]
    RH[count] = RH_mat[j,2+i]
    TMP[count] = TMP_mat[j,2+i]
    count = count + 1
  }
}

dat = data.frame(obs_ind = obs_ind,Date = Date,Hour = Hour , time_ind = time_ind,
                 Station = Station,Region = Region,Lat = Lat,Lon = Lon ,
                 loc_ind = loc_ind,PM10 = PM10,O3 = O3, RH = RH, TMP = TMP)

########## Create a December data

nt_dec = 31*24
dat_dec = matrix(ncol=10,nrow= ns*nt_dec)
colnames(dat_dec) = c("Date","Hour","time_ind","Station","Region","loc_ind","Lon","Lat","PM10","O3")

Date = O3_mat_dec[rep(1:nt_dec,each=ns),"FECHA"]
Hour = O3_mat_dec[rep(1:nt_dec,each=ns),"HORA"]
time_ind = rep((-nt_dec + 1):0,each= ns)
Station = rep(coords[,"Stations"],times=nt_dec)
Region = rep(coords[,"Region"],times=nt_dec)
Lat = rep(coords[,"lat"],times=nt_dec)
Lon = rep(coords[,"lon"],times=nt_dec)
loc_ind = rep(1:ns,times= nt_dec)
PM10 <- O3 <- numeric(nt_dec*ns)
count = 1

for(j in 1:nt_dec){
  for(i in 1:ns){
    PM10[count] = PM10_mat_dec[j,2+i]
    O3[count] = O3_mat_dec[j,2+i]
    count = count + 1
  }
}

dat_dec = data.frame(Date = Date,Hour = Hour , time_ind = time_ind,
                 Station = Station,Region = Region,Lat = Lat,Lon = Lon ,
                 loc_ind = loc_ind,PM10 = PM10,O3 = O3)

#idx_O3 = which(dat_dec$O3 == -99 )
#idx_PM10 = which(dat_dec$PM10 == -99 )

#dat_dec$O3[idx_O3] = mean(dat_dec$O3[-idx_O3] )
#dat_dec$PM10[idx_PM10] = mean(dat_dec$PM10[-idx_PM10])

########## Distance matrix

dd = matrix(0,ns,ns)
for(i in 1:ns){
  dd[i,] = hd(coords$lon[i],coords$lat[i], coords$lon, coords$lat)
}

########## Proximity matrix


W = dd / max(dd)
W = exp(-W)
diag(W) = 0
w_plus = apply(W,1,sum)
Dw = diag(w_plus)
Q = Dw - W

cor(dat$PM10,dat$O3)
cor(dat$PM10,dat$RH)
cor(dat$PM10,dat$TMP)
cor(dat$PM10,dat$TMP * dat$RH)

cor(dat$O3,dat$RH)
cor(dat$O3,dat$TMP)

acf(dat$O3[dat$loc_ind == 1],lag.max = 7*24)
abline(v = 12,lty=3,col="red")
abline(v = 24,lty=3,col="red")

pacf(dat$O3[dat$loc_ind == 1],lag.max = 7*24)
abline(v = 24,lty=3,col="red")

acf(dat$PM10[dat$loc_ind == 1],lag.max = 7*24)
abline(v = 12,lty=3,col="red")
abline(v = 24,lty=3,col="red")

pacf(dat$PM10[dat$loc_ind == 1],lag.max = 7*24)
abline(v = 24,lty=3,col="red")

n_hold  = round(nt*ns*0.1 )

rm(list=setdiff(ls(), c("dat","dat_dec","Q","ns","nt","n_hold","nt_dec")))

