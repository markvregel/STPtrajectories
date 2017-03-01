<img src="logo.png" align="right" width="50%" height="50%"/>

Overview
--------

PACKAGE IS INCOMPLETE AND CONTAINS ERRORS. Package for handling Space-Time Prism(STP) trajectories.
 It contains functions to calculate Potential Path Areas(PPAs), create random
 trajectories and to test for possible encounters by applying the alibi query.
 It also provides functions to visulize the STPs treajectories in 3D.
Installation
------------

``` r

# the the development version from GitHub:
# install.packages("devtools")
devtools::install_github("markvregel/STPtrajectories")
```

Usage
-----


``` r
library(spacetime)
library(sp)
#--------------------------create a STP_Track--------------------------
## create trajectory data
t1 <- as.POSIXct(strptime("01/01/2017 12:00:00", "%m/%d/%Y %H:%M:%S"))
t2 <- as.POSIXct(strptime("01/7/2017 12:00:00", "%m/%d/%Y %H:%M:%S"))
time <- seq(t1,t2,2*60*60)
n <- length(time)
x = cumsum(runif(n) * 8000)
y = smooth(cumsum(runif(n,-0.7,1) * 16000))

crs_NL = CRS("+init=epsg:28992")

points <- SpatialPoints(cbind(x,y),crs_NL)

temp <-18 + cumsum(runif(n,-0.3,0.25))
altitude <- 200 + cumsum(runif(n,-0.75,1)*50)

data <- data.frame(temperature = temp, elevation = altitude)
## create a STP_track
# create class STIDF
stidf1 = STIDF(points, time, data)

# Track-class {trajectories}
my_track1<-Track(stidf1)

# set maximum speed
v1<-10/3.6# speed 10 km/h = 2.777778 m/s

# STP_track class
STP_track1<-STP_Track(my_track1,v1)
# plot
plot(STP_track1,type='p',pch=19,cex=0.8)
# calculate PPA and add to plot
PPA<-calculate_PPA(STP_track1)
plot(PPA,add=TRUE)
```

[Also take look at package trajectories](https://github.com/edzer/trajectories)

Learning to use STPtrajectories
----------------

Getting help
------------
1. Check function help
2. Search in the manual: STPtrajectories.pdf
3. Read the vignettes

