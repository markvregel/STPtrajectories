library(STPtrajectories)
require(spacetime)
## create a STP_track
t1 <- strptime("01/01/2017 00:00:00", "%m/%d/%Y %H:%M:%S")
t2 <- t1+5*60*60
time<-seq(t1,t2,30*60)

x=c(seq(0,25,5),seq(27.5,37.5,2.5))
y=sample(-2:2, 11,replace = T)

n = length(x)
crs_NL = CRS("+init=epsg:28992")

# create class STIDF
stidf1 = STIDF(SpatialPoints(cbind(x,y),crs_NL), time, data.frame(co2 = rnorm(n),O2=rnorm(n)))

# Track-class {trajectories}
my_track1<-Track(stidf1)

# set maximum speed
v1<-getVmaxtrack(my_track1)+0.00015

# STP_track class
STP_track1<-STP_Track(my_track1,v1)

## PPA entire track
#calculate PPA
PPA<-calculate_PPA(STP_track1)

# plot results
plot(STP_track1,type='b')
plot(PPA,add=T)

## PPA only using every second point
# calculate PPA
PPA<-calculate_PPA(STP_track1,points = seq(1,11,2))

# plot results
plot(STP_track1,type='b')
plot(PPA,add=T)
## PPA of a specfic moment in time
# calculate PPA
time <- strptime("01/01/2017 01:15:00", "%m/%d/%Y %H:%M:%S")
PPA<-calculate_PPA(STP_track1,time = time)

# plot results
plot(STP_track1,type='b')
plot(PPA,add=T)

## PPA for a time range
# calculate PPA
timerange1 <- c(t1,strptime("01/01/2017 02:15:00", "%m/%d/%Y %H:%M:%S"))
PPA<-calculate_PPA(STP_track1,time = timerange1)

# plot results
plot(STP_track1,type='b')
plot(PPA,add=T)



