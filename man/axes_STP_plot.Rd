% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/STP_plot.R
\name{axes_STP_plot}
\alias{axes_STP_plot}
\title{axes_STP_plot}
\usage{
axes_STP_plot(minmaxT, z_factor, n_ticks_xy = 3, n_ticks_z = 5,
  expand = 1.1)
}
\arguments{
\item{minmaxT}{a vector of length 2 with two "POSIXct" or"POSIXt" values.
The first and last moment in time of plotted tracks.
Make sure first time is equal to the tracks that starts as first.
Also take into account time_uncetrainty.}

\item{z_factor}{the z facfor used in the plot}

\item{n_ticks_xy}{number of ticks used for the x and y axes}

\item{n_ticks_z}{number of ticks used for the z axes}

\item{expand}{for expanding the bbox. passed to rgl.bbox.
If not all z/time tick are visisble, increase expand value.}
}
\description{
This function adds a bbox with axis to a \link{STP_plot}
}
\examples{
## create a STP_track of two points/ one Space-time Prism(STP)
# time
t1 <- as.POSIXct(strptime("01/01/2017 00:00:00", "\%m/\%d/\%Y \%H:\%M:\%S"))
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
open3d()
z_fac<-2 # relative size of z scale/aspect ratio to spatial scale
# plot STP
STP_plot(STP1,time_interval = 0.5,zfactor = z_fac)


# add axes
axes_STP_plot(time,z_factor = z_fac)
}
\author{
Mark ten Vregelaar
}
