#' apply functions and visualise results
#' BELANGERIJK: verbeter subest functie. dataframe één colom gaat fout en kan zoeizo beter

# clear workspace
rm(list = ls())

#load libraries

library(STPtrajectories)
library(spacetime)
# library(microbenchmark)<--------------



# get bear data
load("d_ata/mating2009.Rdata")
mating2009$LMT_date <- as.POSIXct(mating2009$LMT_date,
                                  format='%d-%m-%Y %H:%M:%S')
bearnames <- sort(unique(mating2009$Name))
# one track
Abborgina <- subset(mating2009, Name==bearnames[1])
Abborgina <- Abborgina[order(Abborgina$LMT_date),]


# put into instance of track class
crs <- CRS("+init=epsg:2400")
keeps <- c("Longitude", "Latitude",'DOP','LMT_date')
Abbo50 <- Abborgina[0:200,]
Abbo50_rd <- Abbo50[keeps]
Abbo50_track_st <- STIDF(SpatialPoints(cbind(Abbo50$Locale_E, Abbo50$Locale_N),crs),Abbo50$LMT_date,Abbo50_rd)
Abbo50_track <- Track(Abbo50_track_st)

# put into instance of STP class

STP1 <- STP_Track(Abbo50_track,rep(0.444444444,length(Abbo50_track)-1))# speed in km/h #1.6 km/h vmax
# only first 50 points
Abbo10_STP_track <- STP1[1:50,'']
#only bewtween time of first 10 points
Abbo10_STP_track <- STP1[1:length(STP1),paste(STP1@endTime[1],STP1@endTime[50],sep="::")]

#STP1[STP1@sp,STP1@time[1:50]]
STPt<-STP1@endTime[19]+15*60

# calc PPA for specific time
PPA <-calculate_PPA(Abbo10_STP_track,STPt)
# plot resulting PPA and trajectory
plot(STP1,type='b')
plot(PPA,add=T,col=rgb(0, 1, 0,0.5))

# calc PPA entire track
STP_moving<-STP1[67:75,'']
# set max speed
STP_moving@connections$vmax<-getVmaxtrack(STP_moving)+0.0005

# PPA for specific time-interval
t1<-STP_moving@endTime[1]+12*60
t2<-STP_moving@endTime[4]+14*60

# Start the clock!
ptm <- proc.time()
result<-calculate_PPA(STP_moving,time = c(t1,t2),x_density = 250,time_interval = 0.5,quadsegs = 12)

# Stop the clock
proc.time() - ptm
plot(result)
## comparing two PPA methods
#source('R/PPA_functions.R')
# Start the clock!
ptm <- proc.time()
result1<-calcPPA_STP_Track(STP_moving)
# Stop the clock
proc.time() - ptm
plot(result1,border='blue')

#source('R/PPA_functions2.R')
# Start the clock!
ptm <- proc.time()
result2<-calculate_PPA(STP_moving,x_density = 1000)
# Stop the clock
proc.time() - ptm
plot(result2,border='red')



# tests:
calculate_PPA(STP_moving,index(STP_moving@time[1]),1:2)
result<-calculate_PPA(STP_moving,c(STP_moving@endTime[1],t2),time_interval = 0.5)
plot(result)

result<-calculate_PPA(STP_moving,c(STP_moving@endTime[1],index(STP_moving@time[1])+7*60),time_interval=0.5)
plot(result)

result<-calculate_PPA(STP_moving,c(STP_moving@endTime[2]-15*60,index(STP_moving@time[3])+15*60),time_interval=0.5)
plot(result)

result<-calculate_PPA(STP_moving,c(STP_moving@endTime[2],index(STP_moving@time[4])),time_interval=0.5)
plot(result)

result<-calculate_PPA(STP_moving,c(STP_moving@endTime[2]-15*59,index(STP_moving@time[2])+15*59),time_interval=0.5)
plot(result)

result<-calculate_PPA(STP_moving,c(index(STP_moving@time[2])-15*60,index(STP_moving@time[2])+15*60),time_interval=0.5)
plot(result)

result<-calculate_PPA(STP_moving,c(STP_moving@endTime[2],index(STP_moving@time[3])),time_interval=0.5)
plot(result)

result<-calculate_PPA(STP_moving,c(STP_moving@endTime[2],index(STP_moving@time[2])+0.5*60),time_interval=0.5)
plot(result,add=T)

result<-calculate_PPA(STP_moving,c(STP_moving@endTime[2]+5*60,index(STP_moving@time[2])+10*60),time_interval=0.5)
plot(result,add=T)

# gives error because time_intercal is too small
result<-calculate_PPA(STP_moving,c(index(STP_moving@time[1]),index(STP_moving@time[1])+0.4*60),time_interval=0.5)
plot(result)
#

result<-calculate_PPA(STP1,c(index(STP1@time[1]),index(STP1@time[99])),time_interval=0.5)
plot(result)
plot(STP1[1:99,''],add=T)
tracks = Tracks(list(tr1 = STP1[1:99,'']))
tracksCollection = TracksCollection(list(tr = tracks))
stplot(tracksCollection, attr = "speed", lwd = 3, scales = list(draw = F),add=T)



## alibi query


# query for two STPs
# first bead
t1=0.0; x1=-1.0; y1=-1.0;
t2=1.0;x2=-1.5; y2=-1.0;
v1=4/60/60;# 4.0
# second bead
t3=0.0; x3=1.0; y3=1.0;
t4=2.0; x4=1.0; y4=1.0;
v2=2.26/60/60



t1 <- as.POSIXct(as.Date("2013-09-30",tz="CET"))
t2 <- t1+t2*60*60;t3 <- t1;t4 <- t3 +t4*60*60

xa = c(x1,x2);ya = c(y1,y2)
xb = c(x3,x4);yb = c(y3,y4)

n = length(xa)
crs_NL = CRS("+init=epsg:28992")

# create class STIDF
stidf1 = STIDF(SpatialPoints(cbind(xa,ya),crs_NL), c(t1,t2), data.frame(co2 = rnorm(n)))
stidf2 = STIDF(SpatialPoints(cbind(xb,yb),crs_NL), c(t3,t4), data.frame(co2 = rnorm(n)))


# Track-class {trajectories}
my_track1<-Track(stidf1)
my_track2<-Track(stidf2)

# STP_track class
STP_1<-STP_Track(my_track1,v1)
STP_2<-STP_Track(my_track2,v2)

alibi_STP(STP_1,STP_2)

STP_plot(STP_1,0.4,225/60/60)
STP_plot(STP_2,0.4,225/60/60,'blue',STP_1@endTime[1])
axes_STP_plot(STP_2@endTime,225/60/60,3,5)
# query for entire trajcetories

# Read August 2007 data
data2007 <- read.table("d_ata/August2007.txt", header = T, sep=",")
data2007$LMT_date <- as.POSIXct(data2007$LMT_date,
                                format='%d-%m-%Y %H:%M:%S')
plot(data2007$Locale_E, data2007$Locale_N, asp=1)
# delete obvious outlier caused by GPS error
data2007 <- subset(data2007, Locale_E > 1400000)
plot(data2007$Locale_E, data2007$Locale_N, asp=1)

# SUBSETS and ordering
Koski2007 <- subset(data2007, PubName == "Koski (2310)")
Grivla2007 <- subset(data2007, PubName == "Grivla (2911)")
Koski2007 <-Koski2007[order(Koski2007$LMT_date),] # order on data_time
Grivla2007 <- Grivla2007[order(Grivla2007$LMT_date),]
# 2D Plot trajectories
plot(Koski2007$Locale_E, Koski2007$Locale_N, type="l",
     main = "Track Koski, August 2007") # 2D single bear

# plot a space-time cube
tstart <- min(min(Koski2007$LMT_date),min(Grivla2007$LMT_date))
Koski2007$tspan <- difftime(Koski2007$LMT_date, tstart, units="hours")
Grivla2007$tspan <- difftime(Grivla2007$LMT_date, tstart,
                             units="hours")
xscale <- c(min(Grivla2007$Locale_E,Koski2007$Locale_E),
            max(Grivla2007$Locale_E,Koski2007$Locale_E))
yscale <- c(min(Grivla2007$Locale_N,Koski2007$Locale_N),
            max(Grivla2007$Locale_N,Koski2007$Locale_N))
zscale <- c(0,max(Grivla2007$tspan,Koski2007$tspan))
open3d(windowRect=c(100,100,700,700))
plot3d(Koski2007$Locale_E, Koski2007$Locale_N, Koski2007$tspan,
       type="l", col="red", lwd=1.5,
       xlab="East", ylab="North", zlab="time", xlim=xscale,
       ylim=yscale,zlim=zscale)
lines3d(Grivla2007$Locale_E, Grivla2007$Locale_N, Grivla2007$tspan,
        col="blue", lwd= 1.5)





Koski_rd <- Koski2007[keeps]
Koski_track_st <- STIDF(SpatialPoints(cbind(Koski2007$Locale_E, Koski2007$Locale_N),crs),Koski2007$LMT_date,Koski_rd)
Koski_track <- Track(Koski_track_st)

# put into instance of STP class
vmax_Koski<-getVmaxtrack(Koski_track)+0.277777778
Koski_STP_track <- STP_Track(Koski_track,vmax_Koski)

Grivla_rd <- Grivla2007[keeps]
Grivla_track_st <- STIDF(SpatialPoints(cbind(Grivla2007$Locale_E, Grivla2007$Locale_N),crs),Grivla2007$LMT_date,Grivla_rd)
Grivla_track <- Track(Grivla_track_st)

# put into instance of STP class
vmax_Grivla<-getVmaxtrack(Grivla_track)+0.277777778
Grivla_STP_track <- STP_Track(Grivla_track,vmax_Grivla)

points<-c()
# remove points that are within two minutes
for (i in 2:length(Koski_STP_track)){
  if(difftime(Koski_STP_track@endTime[i],
              Koski_STP_track@endTime[i-1],units = 'secs')<=120){

    points<-c(points,i)
  }
}

Koski_track<-Koski_STP_track[Koski_STP_track@sp[(points*-1)],'']
vmax_Koski<-getVmaxtrack(Koski_track)+0.277777778
Koski_STP_track <- STP_Track(Koski_track,vmax_Koski)
points<-c()
# remove points that are within two minutes
for (i in 2:length(Grivla_STP_track)){
  if(difftime(Grivla_STP_track@endTime[i],
              Grivla_STP_track@endTime[i-1],units = 'secs')<=120){

    points<-c(points,i)
  }
}

Grivla_track<-Grivla_STP_track[Grivla_STP_track@sp[(points*-1)],'']
vmax_Grivla<-getVmaxtrack(Grivla_track)+0.277777778# 1km/h
Grivla_STP_track <- STP_Track(Grivla_track,vmax_Grivla)


#subset trajectories
Koski_STP_track<-Koski_STP_track[125:250,'']
Grivla_STP_track<-Grivla_STP_track[125:550,'']


## visulise in 3D

#STP_moving<-subset_STP_Track(STP_moving,1:6)
n<-length(STP_moving)
zf<-60
STP_moving@connections$vmax<-0.60
y<-(bbox(STP_moving)[2,1:2])
y[2]<-y[2]+500
y[1]<-y[1]-500


open3d()
stcube(STP_moving,ylim=y,showMap=F,type='p',lwd=5,zlab='',xlab='',ylab='',size=5,col='blue')# ----- plot with basemap

#STP_Track.plot(STP_moving,0.35,zf)

minmaxTime<-c(STP_moving@endTime[1],STP_moving@endTime[n])

tdif<-as.numeric(difftime(minmaxTime[2],minmaxTime[1],units = 'mins'))
tickval<-seq(0,tdif*zf,tdif*zf/5)

timesval<-seq(minmaxTime[1],minmaxTime[2],(tdif*60)/5)

rgl.bbox(yat=0,xat = 0,zat = tickval, zlab = timesval,color = c("black", "black"), emission = "#252526",
         specular = "#636363", shininess = 5, alpha = 0.8,expand=1.02,marklen = 5)

axes3d(c('x--'),nticks = 3,texmipmap=T)
axis3d(c('y-+'),at = c(6811500,6812500),labels=c(6811500,6812500),ntick=2,line = 2)
title3d(main='Small Track', zlab='time',xlab = 'x(meters)',color='black',size=5,line=2.5)
mtext3d('y(meters)', 'y-+',line =0.5)
bg3d('white')
STP_plot(STP_moving,1.2,zf)

axes3d(c('z++'),at = tickval, labels = timesval)
# crazy
# STP_crazy<-STP_moving
# STP_crazy@endTime<-STP_crazy@endTime-15*60
# STP_Track.plot(STP_crazy,1.2,zf,'blue',STP_moving@endTime[[1]])




########### plot STP_tracks based on time interval
zf<-17
Koski_STP_track <- STP_Track(Koski_track,vmax_Koski)
Grivla_STP_track <- STP_Track(Grivla_track,vmax_Grivla)
t1<-strptime("2007-08-21 02:00:00", "%Y-%m-%d %H:%M:%S")
t2<-strptime("2007-08-21 22:00:00", "%Y-%m-%d %H:%M:%S")


STP_track1b<-Koski_STP_track[1:9999,paste(t1,t2,sep = '::')]
STP_track2b<-Grivla_STP_track[1:9999,paste(t1,t2,sep = '::')]


addspeed<-3#0.3
STP_track1b@connections$vmax<-getVmaxtrack(STP_track1b)+addspeed
STP_track2b@connections$vmax<-getVmaxtrack(STP_track2b)+addspeed


#look at threading-----------------------------------------
open3d()
STP_plot(STP_track1b,2,zf)
STP_plot(STP_track2b,2,zf,'blue',STP_track1b@endTime[[1]])
minmaxTime<-c(STP_track1b@endTime[1],tail(STP_track1b@endTime,1))
axes_STP_plot(minmaxTime,zf,7)

alibiQ_STP_Track(STP_track1b,STP_track2b)
##### plot something two other tracks,little movement
STP_track1<-Koski_STP_track[125:250,'']
STP_track2<-Grivla_STP_track[125:250,'']

STP_track1a<-STP_track1[40:45,'']
STP_track2a<-STP_track2[64:71,'']

zf=200
open3d()
STP_plot(STP_track1a,0.8,zf)
STP_plot(STP_track2a,0.8,zf,'blue',st=STP_track1a@endTime[1])

minmaxTime<-c(STP_track1a@endTime[1],tail(STP_track1a@endTime,1))
axes_STP_plot(minmaxTime,zf,5)



