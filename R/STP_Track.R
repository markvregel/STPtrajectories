#' @title STP_Track class
#' @description A class to represent Space-Time Prism(STP) trajecories.
#' These are trajectories with a maximum speed for each segment.
#' The maximum speed is added to the connections slot
#' @param Track object of class  \link{Track}
#' @inheritSection trajectories::Track Slots of class "Track"
#' @examples
#'library(spacetime)
#'library(sp)
#'#--------------------------create a STP_Track--------------------------
#'#------------------------------example 1------------------------------
#'## create trajectory data
#'t1 <- as.POSIXct(strptime("01/01/2017 12:00:00", "%m/%d/%Y %H:%M:%S"))
#'t2 <- as.POSIXct(strptime("01/7/2017 12:00:00", "%m/%d/%Y %H:%M:%S"))
#'time <- seq(t1,t2,2*60*60)
#'n <- length(time)
#'x = cumsum(runif(n) * 8000)
#'y = smooth(cumsum(runif(n,-0.7,1) * 16000))
#'
#'crs_NL = CRS("+init=epsg:28992")
#'
#'points <- SpatialPoints(cbind(x,y),crs_NL)
#'
#'temp <-18 + cumsum(runif(n,-0.3,0.25))
#'altitude <- 200 + cumsum(runif(n,-0.75,1)*50)
#'
#'data <- data.frame(temperature = temp, elevation = altitude)
#'## create a STP_track
#'# create class STIDF
#'stidf1 = STIDF(points, time, data)
#'
#'# Track-class {trajectories}
#'my_track1<-Track(stidf1)
#'
#'# set maximum speed
#'v1<-10/3.6# speed 10 km/h = 2.777778 m/s
#'
#'# STP_track class
#'STP_track1<-STP_Track(my_track1,v1)
#'# plot
#'plot(STP_track1,type='p',pch=19,cex=0.8)
#'# calculate PPA and add to plot
#'PPA<-calculate_PPA(STP_track1)
#'plot(PPA,add=TRUE)
#'
#'#------------------------------example 2------------------------------
#'## vmax depends on elevation
#'# assuming that max speed is lower as result of the thinner air
#'vmax<- getVmaxtrack(STP_track1) + (max(STP_track1@data$elevation[1:n-1])-STP_track1@data$elevation[1:n-1])*0.01
#'STP_track1@connections$vmax<-vmax
#'# calculate PPA
#'PPA<-calculate_PPA(STP_track1)
#'# create tracksCollection and plot
#'tracks = Tracks(list(tr1 = STP_track1))
#'tracksCollection = TracksCollection(list(tr = tracks))
#'stplot(tracksCollection, attr = "elevation", lwd = 3, scales = list(draw =TRUE),
#'       sp.layout=PPA,xlim = PPA@bbox[1,],ylim = PPA@bbox[2,],main= "Track with PPA\n vmax depends on altitude",
#'       sub='colour is altitude in meters',xlab='x',ylab='y')
#'#------------------------------example 3------------------------------
#'## vmax depends on the distance to get to next point.
#'# Thus on the distance that needs to be covered in the avialable time
#'# Assuming that if two points are closer together the max speed is lower
#'STP_track1@connections$vmax<-getMinimalSpeed(STP_track1)*1.5
#'# calculate PPA
#'PPA<-calculate_PPA(STP_track1)
#'# create tracksCollection and plot
#'plot(PPA,add=TRUE)
#'tracks = Tracks(list(tr1 = STP_track1))
#'tracksCollection = TracksCollection(list(tr = tracks))
#'stplot(tracksCollection, attr = "vmax", lwd = 3, scales = list(draw = TRUE),
#'       sp.layout=PPA,xlim = PPA@bbox[1,],ylim = PPA@bbox[2,],
#'       main= "Track with PPA\n vmax depends on the distance and time budget between consecutive points",
#'       sub='colour is vmax in m/s',xlab='x',ylab='y',cex.main = 0.75)
#'
#'
#'#--------------------------subset a STP_Track--------------------------
#'# make sure vmax is high enough
#'STP_track1@connections$vmax<-getVmaxtrack(STP_track1)+1
#'#------------------------------example 1------------------------------
#'## subset based on space-time points
#'# get first 10 points
#'STP_track1_10<-STP_track1[1:10,'']
#'# only keep every second point
#'STP_track1_depleted<-STP_track1[seq(1,n,2),'']
#'plot(STP_track1,type='b')
#'plot(STP_track1_depleted,type='b')
#'#------------------------------example 2------------------------------
#'## subset based on time
#'STP_track1_a<-STP_track1[1:n,'2017-01-01 12:00:00 CET::2017-01-02 16:00:00 CET']
#'# only keep every second point within the time interval
#'STP_track1_b<-STP_track1[seq(1,n,2),'2017-01-01 12:00:00 CET::2017-01-02 16:00:00 CET']
#'# all points in the night(between 22:00 and 08:00)
#'STP_track1_c<-STP_track1[1:n,"T22:00/T08:00"] # see package xts for more handy subsetting tricks
STP_Track<-setClass(
  # Class for storing Space-time Prism trajecotries. These are trajecories with maximum speeds.
  "STP_Track",
  contains = "Track",

  validity = function(object) {



    speedCheck<-object@connections$vmax>=getMinimalSpeed(object)


    if((FALSE %in% speedCheck)) {
      return(paste0("The vmax for connection ",which(speedCheck %in% F),
                    " is smaller than minimal speed required to reach the next point"))
    }
    # if not projected not a valid STP track
    if(is.projected(object)){
      return(TRUE)
    }
    else{
      return('Track is not projected, unable to create spac-time prisms.')
    }
  }
)
#'@export
STP_Track = function(track,vmax) {
  # vmax distance unit projection/seconds
  # degrees are converted to Great Circle distance in Trmeters and thus m/s
  track@connections$vmax<- vmax
  validObject(track)
  new("STP_Track", track)
}

#
subs.STP_Track <- function(x, i, j, ..., drop = TRUE) {
  # Provide selection methods.
  # i selection of record index (spatial/temporal/spatio-temporal entities)
  # j selection of temporal entities (see syntax in package xts)
  track<-Track(as(x, "STIDF")[i, j, ..., drop = drop])

  vmax<-x[['vmax']][x@endTime %in% track@endTime]
  if (max(diff(which(x@endTime %in% track@endTime,TRUE)))>1){
    warning('Some intermediate space-time points will be removed. If vmax values differ along the STP_track, the vmax values might change.
            The maxmimum speed of a new line segments will depend on the first point of the segment')
  }
  return(STP_Track(track,head(vmax,-1)))

}

setMethod("[", "STP_Track", subs.STP_Track)




