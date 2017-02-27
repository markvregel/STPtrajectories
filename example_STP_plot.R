library(STPtrajectories)
library(spacetime)
## create 2 STP_tracks
# time
t1 <- strptime("01/01/2017 00:00:00", "%m/%d/%Y %H:%M:%S")
t2 <- t1+5*60*60 # 5 hours after t1
time1<-seq(t1,t2,30*60)
time2<-time1+0.25*60*60
# spatial coordinates
x1=c(seq(0,25,5),seq(27.5,37.5,2.5))
y1=sample(-2:2, 11,replace = T)
x2=c(seq(0,25,5),seq(27.5,37.5,2.5))
y2=sample(-2:2, 11,replace = T)

n = length(x1)
crs_NL = CRS("+init=epsg:28992")

# create class STIDF
stidf1 = STIDF(SpatialPoints(cbind(x1,y1),crs_NL), time1, data.frame(co2 = rnorm(n),O2=rnorm(n)))
stidf2 = STIDF(SpatialPoints(cbind(x2,y2),crs_NL), time2, data.frame(co2 = rnorm(n),O2=rnorm(n)))

# Track-class {trajectories}
my_track1<-Track(stidf1)
my_track2<-Track(stidf2)
# set maximum speed
v1<-getVmaxtrack(my_track1)+0.00015
v2<-getVmaxtrack(my_track2)+0.00030
# STP_track class
STP_track1<-STP_Track(my_track1,v1)
STP_track2<-STP_Track(my_track2,v2)


## 3D STP plot of STP_tracks
z_fac<-0.2 # relative size of z scale/aspect ratio to spatial scale
# plot STPS first STP_track
STP_plot(STP_track1,time_interval = 1,z_fac)
# plot STPS second STP_track
STP_plot(STP_track2,time_interval = 1,z_fac,'blue',st = STP_track1@endTime[1])# provide st for correct starting location first STP

# calculate first and last moment in time
min_max_Time<-c(STP_track1@endTime[1],STP_track2@endTime[length(STP_track2)])
# add axes
axes_STP_plot(min_max_Time,z_factor = z_fac)

# add title and change background colour
library(rgl)
title3d(main = '2 randomly generated STP tracks')
bg3d('lightblue')

# example axes_STP_plot
## create a STP_track of two points/ one Space-time Prism(STP)
# time
t1 <- as.POSIXct(strptime("01/01/2017 00:00:00", "%m/%d/%Y %H:%M:%S"))
t2 <- t1+0.5*60*60 # 30 after t1
time<-c(t1,t2)
# spatial coordinates
x=c(0,2);y=c(1,3)

n = length(x)
crs_NL = CRS("+init=epsg:28992")

# create class STIDF
stidf1 = STIDF(SpatialPoints(cbind(x,y),crs_NL), time, data.frame(co2 = rnorm(n),O2=rnorm(n)))

# Track-class {trajectories}
my_track1<-Track(stidf1)

# set maximum speed
v1<-getVmaxtrack(my_track1)+0.01
# STP_track class
STP1<-STP_Track(my_track1,v1)

## 3D STP plot of STP

z_fac<-1.5 # relative size of z scale/aspect ratio to spatial scale
# plot STP
STP_plot(STP1,time_interval = 0.5,z_fac)


# add axes
axes_STP_plot(time,z_factor = z_fac)
