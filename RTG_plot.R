#' @title RTG
#' @description Random Trajecory generator
#'
#' @param n_points number of points to add between two consecutive space-time points
#' @param STP_track the  \link{STP_Track} to which the randomly generated space-time points are added
#' @param max_time_interval  The max_time_interval determines between which space-time points random points are added.
#' If the time difference between two points is bigger than max_time_interval,
#' n_points number of points will be added between the two points.
#' If no value is provided new point(s) will be added between all consecutive space-time points
#' @param quadsegs Passed to buffer. Number of line segments to use to approximate a quarter circle.
#' Only used where paramter time_interval is relavant
#' @param iter number of times to try to place sample points in the PPA before giving up and returning NULLter (default = 4) -
#' this may occur when trying to hit a small and awkwardly shaped polygon in a large bounding box with a small number
#' of points.
#' @importFrom rgeos gBuffer
#' @importFrom spacetime STIDF
#' @return a \link{STP_Track}
#' @export
#' @author Mark ten Vregelaar

RTG_plot<-function(STP_track,n_points=1,max_time_interval=NULL,quadsegs=12,iter=4){
  plot(STP_track,type='b')
  # get the segments that exceed the max_time_interval
  if (!is.null(max_time_interval)){
    time_steps<-diff(STP_track@endTime)
    exceed_seg<-which(time_steps>max_time_interval)
    # if all segments within max_time_interval stop
    if(length(exceed_seg)==0){
      stop("All the time difference between consecutive space-time points are smaller than max_time_interval")
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
    startpoint <- SpatialPointsDataFrame(STP@sp[1,],data=data.frame(time=t1))
    endpoint <- SpatialPointsDataFrame(STP@sp[2,],data=data.frame(time=t2))
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
      plot(PPA,add=T,border='red')
      # randomly select point in PPA
      npoint<-spsample(PPA,1,type = 'random',iter=iter)
      npoint<-SpatialPointsDataFrame(npoint,data=data.frame(time=t))
      plot(npoint,add=T)
      Sys.sleep(1)
      # add point to the SpatialPointsDataFrame with the two orginal points and other random points
      points<-spRbind(points,npoint)
    }
    # add the newly created point(s) to all_points
    all_points<-rbind(all_points,points[3:n_points,])
  }
  # create class STIDF
  random_STIDF = STIDF(all_points, all_points$time, data.frame(STP@data))
  # Track-class {trajectories}
  random_track<-Track(random_STIDF)
  # STP_track class
  random_STP<-STP_Track(random_track,STP@connections$vmax)

}
