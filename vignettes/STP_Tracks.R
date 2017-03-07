## ----include=FALSE-------------------------------------------------------
library(STPtrajectories)

## ----setup---------------------------------------------------------------
library(sp)
library(spacetime)
library(knitr)
library(rgl)
knit_hooks$set(webgl = hook_webgl)

## ------------------------------------------------------------------------
crs_NL = CRS("+init=epsg:28992")
t1 <- as.POSIXct(strptime("01/01/2017 12:00:00", "%m/%d/%Y %H:%M:%S"))
t2 <- as.POSIXct(strptime("01/01/2017 15:00:00", "%m/%d/%Y %H:%M:%S"))
time<-c(t1,t2)


# Spatial coordinates
x_bear1=c(7500,10);y_bear1=c(10,6000)

x_bear2=c(20,8000);y_bear2=c(30,7000)

n = length(x_bear1)

# create class STIDF {spacetime}
stidf_bear1 = STIDF(SpatialPoints(cbind(x_bear1,y_bear1),crs_NL), time, data.frame(co2 = rnorm(n),O2=rnorm(n)))
stidf_bear2 = STIDF(SpatialPoints(cbind(x_bear2,y_bear2),crs_NL), time, data.frame(co2 = rnorm(n),O2=rnorm(n)))

# Track-class {trajectories}
track_bear1<-Track(stidf_bear1)
track_bear2<-Track(stidf_bear2)
# Set maximum speed 
# max speed for a period of 3 hours
v_bear1<-5/3.6
v_bear2<-6/3.6
# STP_track class
STP_bear1<-STP_Track(track_bear1,v_bear1)
STP_bear2<-STP_Track(track_bear2,v_bear2)
# plot the points
plot(STP_bear1,type='p',col='red',pch=16,cex=2,xlim=c(0,9000),ylim=c(0,8000),xlab='x (meters)',
     ylab='y (meters)', main= 'Start and end locations of bear trajecories')
plot(STP_bear2,type='p',col='blue',pch=16,cex=2,add=T)
legend('topright',c('Bear 1','Bear 2'),pch = 16,col=c('red','blue'))
text(c(x_bear1,x_bear2), c(y_bear1,y_bear2), labels=c('start','end'), cex= 1.2,adj= c(0.5,-1))

## ------------------------------------------------------------------------
## Create trajecotires using the RTG
# set seed to create same trajecory
set.seed(10)
STP_track_bear1<-RTG(STP_bear1,n_points = 15)
set.seed(2)
STP_track_bear2<-RTG(STP_bear2,n_points = 15)

# plot results
plot(STP_track_bear1,type='b',col='red',pch=16,cex=0.8,xlim=c(-1000,12000),ylim=c(-500,9000),xlab='x (meters)',ylab='y (meters)', main= 'The two bear trajecories')
plot(STP_track_bear2,type='b',col='blue',pch=16,cex=0.8,add=T)
legend('topright',c('Bear 1','Bear 2'),pch = 16,lty =1, col=c('red','blue'))


## ------------------------------------------------------------------------
# set the new maximum speed. different for every segment
vmax_bear1<-STP_track_bear1@connections$speed*1.5
vmax_bear2<-STP_track_bear2@connections$speed*1.5

STP_track_bear1@connections$vmax<-vmax_bear1
STP_track_bear2@connections$vmax<-vmax_bear2

# calculate Potential Path Area (PPA)
PPA_bear1 <- calculate_PPA(STP_track_bear1)
PPA_bear2 <- calculate_PPA(STP_track_bear2)


# plot results
plot(STP_track_bear1,type='b',col='red',pch=16,cex=0.8,xlim=c(-1000,12000),ylim=c(-500,9000),xlab='x (meters)',ylab='y (meters)', main= 'The two bear trajecories with PPA')
plot(STP_track_bear2,type='b',col='blue',pch=16,cex=0.8,add=T)
legend('topright',c('Bear 1','Bear 2'),pch = 16,lty =1, col=c('red','blue'))


plot(PPA_bear1,add=T)
plot(PPA_bear2,add=T)


## ----alibi_query---------------------------------------------------------
alibi_query(STP_track_bear1,STP_track_bear2)# not always correct. package still in development


## ----STP_plot, webgl=TRUE------------------------------------------------
# zfac<- 50 # aspect ration between sptatial axes and time axis
# t_int <- 0.8 # determines how many PPAs are used to visualise STPs.
# 
# STP_plot(STP_track_bear1,time_interval = t_int,zfactor = zfac)
# STP_plot(STP_track_bear2,time_interval = t_int,zfactor = zfac,st = STP_track_bear1@endTime[1],col = 'blue')
# 
# # add axes
# #axes_STP_plot(time,z_factor = zfac) function to add axes but not suitable for Rmarkdown
# title3d(main = "3D Visualisation of STP_tracks",xlab='x',ylab='y',cex=1.3)
# bg3d('lightblue')
# # data time axis
# tdif<-as.numeric(difftime(time[2],time[1],units = 'mins'))
# tickval<-seq(0,tdif*zfac,length.out = 5)
# timesval<-seq(time[1],time[2],length.out = 5)
# # add axes
# axes3d(c('x','y'),xlab='x')
# axis3d('z',at=tickval,labels = timesval,cex=0.8)
# box3d()

