library(spacetime)
## Create a random trajecory based on a begin and end point
## Create trajectory with only two points
# Time
t1 <- as.POSIXct(strptime("01/01/2017 00:00:00", "%m/%d/%Y %H:%M:%S"))
t2 <- t1+0.5*60*60 # 2 hours after t1
time<-c(t1,t2)
# Spatial coordinates
x=c(5,10);y=c(10,20)
n = length(x)
crs_NL = CRS("+init=epsg:28992")

# create class STIDF
stidf1 = STIDF(SpatialPoints(cbind(x,y),crs_NL), time, data.frame(co2 = rnorm(n),O2=rnorm(n)))

# Track-class {trajectories}
track1<-Track(stidf1)

# Set maximum speed
v1<-getVmaxtrack(track1)+0.001
# STP_track class
STP1<-STP_Track(track1,v1)
plot(STP1,type='p',col='red',pch=16,cex=2)

# Create a random trajectory between the two points
random_STP_track<-RTG(STP1,n_points = 10)
plot(random_STP_track,type='b',add=T)

## Add points to a trajectory with multple points

## Create a STP_track
np <-6 # Number of points orignal track
t1 <- as.POSIXct(strptime("01/01/2017 00:00:00", "%m/%d/%Y %H:%M:%S"))
random1<-cumsum(sample((0.5*60):(2.8*60*60),np))
time<-t1+random1

x=random1/2
y=seq(1,100,length.out = np)

n = length(x)
crs_NL = CRS("+init=epsg:28992")

# Create class STIDF
stidf2 = STIDF(SpatialPoints(cbind(x,y),crs_NL), time, data.frame(co2 = rnorm(n),O2=rnorm(n)))

# Track-class {trajectories}
track2<-Track(stidf2)

# Set maximum speed
v1<-getVmaxtrack(track2)+0.1

# STP_track class
STP_track2<-STP_Track(track2,v1)
# STP_track2 is track with different time intervals between the space-time points.
# The distance between two points increases with the time interval
plot(STP_track2,type='p',col='red',pch=16,cex=2)

## Fill blank spot of trajecotries in two steps
# Add 2 random points in between two sapce-time points that more than 90 minutes apart
filled_track1 <-RTG(STP_track2,n_points = 2,max_time_interval = 120)
plot(filled_track1,type='p',pch=16,add=T,col='blue')

# Add 1 random point in between two sapce-time points that more than 45 minutes apart
filled_track2 <-RTG(filled_track1,n_points = 1,max_time_interval = 60)
plot(filled_track2,type='b',add=T,cex=0.7)



