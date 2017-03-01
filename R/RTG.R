#' @title RTG
#' @description This random trajectory generator (RTG) works as described in Technitis et al.(2015).
#' It creates a trajectory based on a the space-time prism concept,
#' it randomly adds a user defined number of points between the points of a trajectory.
#' The new point is randomly placed in the PPA of the corresponding point in time.
#' The added points are evenly divided over time and are always within the space-time prism.
#'
#' @param STP_track the  \link{STP_Track} to which the randomly generated space-time points are added
#' @param max_time_interval  The max_time_interval determines between which space-time points the random points are added.
#' If the time difference between two points is bigger than max_time_interval,
#' @param n_points number of points will be added between the two points.
#' If no value is provided new point(s) will be added between all consecutive space-time points
#' @param quadsegs Passed to buffer. Number of line segments to use to approximate a quarter circle.
#' Only used where paramter time_interval is relavant
#' @param iter number of times to try to place sample points in the PPA before giving up and returning NULLter (default = 4) -
#' this may occur when trying to hit a small and awkwardly shaped polygon in a large bounding box with a small number
#' of points.
#' @importFrom spacetime STIDF
#' @importFrom maptools spRbind
#' @importFrom rgeos gBuffer
#' @return a \link{STP_Track} with the newly added random space-time points. Slot data has NAs for the new points.
#' Vmax values for new connections are equal to the vmax values of the original connections.
#' @export
#' @author Mark ten Vregelaar
#' @references - Technitis, G., Othman, W., Safi, K., & Weibel, R. (2015).
#' From A to B, randomly: A point-to-point random trajectory generator for animal movement.
#' International Journal of Geographical Information Science, 29(6), 912-934.
#' \url{http://www.tandfonline.com/doi/abs/10.1080/13658816.2014.999682}
#' @examples
#'library(spacetime)
#'library(sp)
#'#------------------------------example 1------------------------------
#'## Create a random trajecory based on a begin and end point
#'## Create trajectory with only two points
#'# Time
#'t1 <- as.POSIXct(strptime("01/01/2017 00:00:00", "%m/%d/%Y %H:%M:%S"))
#'t2 <- t1+0.5*60*60 # 2 hours after t1
#'time<-c(t1,t2)
#'# Spatial coordinates
#'x=c(5,10);y=c(10,20)
#'n = length(x)
#'crs_NL = CRS("+init=epsg:28992")
#'
#'# create class STIDF
#'stidf1 = STIDF(SpatialPoints(cbind(x,y),crs_NL), time, data.frame(co2 = rnorm(n),O2=rnorm(n)))
#'
#'# Track-class {trajectories}
#'track1<-Track(stidf1)
#'
#'# Set maximum speed
#'v1<-getVmaxtrack(track1)+0.001
#'# STP_track class
#'STP1<-STP_Track(track1,v1)
#'plot(STP1,type='p',col='red',pch=16,cex=2)
#'
#'# Create a random trajectory between the two points
#'random_STP_track<-RTG(STP1,n_points = 10)
#'plot(random_STP_track,type='b',add=TRUE)
#'
#'#------------------------------example 2------------------------------
#'## Add points to a trajectory with multple points
#'## Create a STP_track
#'np <-6 # Number of points orignal track
#'t1 <- as.POSIXct(strptime("01/01/2017 00:00:00", "%m/%d/%Y %H:%M:%S"))
#'random1<-cumsum(sample((0.5*60):(2.8*60*60),np))
#'time<-t1+random1
#'
#'x=random1/2
#'y=seq(1,100,length.out = np)
#'
#'n = length(x)
#'crs_NL = CRS("+init=epsg:28992")
#'
#'# Create class STIDF
#'stidf2 = STIDF(SpatialPoints(cbind(x,y),crs_NL), time, data.frame(co2 = rnorm(n),O2=rnorm(n)))
#'
#'# Track-class {trajectories}
#'track2<-Track(stidf2)
#'
#'# Set maximum speed
#'v1<-getVmaxtrack(track2)+0.1
#'
#'# STP_track class
#'STP_track2<-STP_Track(track2,v1)
#'# STP_track2 is track with different time intervals between the space-time points.
#'# The distance between two points increases with the time interval
#'plot(STP_track2,type='p',col='red',pch=16,cex=2)
#'
#'## Fill blank spot of trajecotries in two steps
#'# Add 2 random points in between two sapce-time points that more than 90 minutes apart
#'filled_track1 <-RTG(STP_track2,n_points = 2,max_time_interval = 120)
#'plot(filled_track1,type='p',pch=16,add=TRUE,col='blue')
#'
#'# Add 1 random point in between two sapce-time points that more than 45 minutes apart
#'filled_track2 <-RTG(filled_track1,n_points = 1,max_time_interval = 60)
#'plot(filled_track2,type='b',add=TRUE,cex=0.7)
RTG<-function(STP_track,n_points=1,max_time_interval=NULL,quadsegs=12,iter=4){
  # DATA GOES LOST??!!!!
  # get the segments that exceed the max_time_interval
  if (!is.null(max_time_interval)){
    time_steps<-diff(STP_track@endTime)
    units(time_steps)<-'mins'
    exceed_seg<-which(time_steps>max_time_interval)

    # if all segments within max_time_interval stop
    if(length(exceed_seg)==0){
      warning("All the time difference between consecutive space-time points are smaller than max_time_interval.
               No new points added")
      return(STP_track)
    }
  }else{
    # if no max_time_interval is provided add add point to all segments
    exceed_seg<-1:(length(STP_track)-1)
  }
  # create SpatialPointsDataFrame with all existing points
  all_points<-SpatialPointsDataFrame(STP_track@sp,data=data.frame(time=STP_track@endTime))

  n_points <- n_points +2  # add two for already existing points
  # loop through segments to add points
  for (i in exceed_seg){
    # get STP that corresponds with segment
    STP<-STP_track[i:(i+1),'']
    # obtain time and speed
    t1 <- STP@endTime[[1]]
    t2 <- STP@endTime[[2]]
    v <- STP@connections$vmax

    # Times for which PPA needs to be calculated in random order to avoid drifting(see Technitis et al.)
    PPA_times<-seq(t1,t2,length.out = n_points)
    shuffled_times<-sample(PPA_times[2:(length(PPA_times)-1)])

    # initialize variables for loop through shuffled_times
    startpoint <- SpatialPointsDataFrame(STP@sp[1,],data=data.frame(time=t1),match.ID = FALSE)
    endpoint <- SpatialPointsDataFrame(STP@sp[2,],data=data.frame(time=t2),match.ID = FALSE)
    points<-rbind(startpoint,endpoint)#SpatialPointsDataFrame with orignal space-time points

    for (j in 1:length(shuffled_times)){
      t<-shuffled_times[j]# the time for which a new point needs to be added

      # get the two surrounding space-time points for the newly to add point
      t1<-max(points$time[points$time<=t])
      t2<-min(points$time[points$time>=t])
      startpoint<-points[points$time==t1,]
      endpoint <- points[points$time==t2,]

      # get time difference between the two surrounding points and t
      t_dif1 = abs(difftime(t1,t,units = 'secs'))
      t_dif2 = abs(difftime(t2,t,units = 'secs'))

      # calculate the maximum travel dictance starting form the two points
      s1 = v*as.numeric(t_dif1)
      s2 = v*as.numeric(t_dif2)

      # calculate the PPA based on itersections of the two buffers
      buffer1<-gBuffer(startpoint,width=s1,quadsegs=quadsegs)
      buffer2<-gBuffer(endpoint,width=s2,quadsegs=quadsegs)
      PPA<-gIntersection(buffer1,buffer2)

      # randomly select point in PPA
      npoint<-spsample(PPA,1,type = 'random',iter=iter)
      npoint<-SpatialPointsDataFrame(npoint,data=data.frame(time=t))

      # add point to the SpatialPointsDataFrame with the two orginal points and other random points
      points<-spRbind(points,npoint)# problems with rbind
    }
    # add the newly created point(s) to all_points
    all_points<-rbind(all_points,points[points$time %in% shuffled_times,])# points[3:n_points,]

  }
  # create class STIDF
  random_STIDF = STIDF(all_points, all_points$time, data.frame(STP@data))
  # return original data
  random_STIDF@data[random_STIDF@endTime %in% STP_track@endTime,] = STP_track@data
  # Track-class {trajectories}
  random_track<-Track(random_STIDF)
  # STP_track class
  random_STP<-STP_Track(random_track,STP@connections$vmax)

}

