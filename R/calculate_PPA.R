## Functions for calculating the Potetial Path Area
#' @title calculate_PPA
#' @description Function for calculating the Potetial Path Area(PPA) of a STP_track.
#' This function can calculate the PPA for the entire trajectory, a specfic moment in time or a time range.
#' @param STP_track The STP_track for which the PPA needs to be calculated
#' @param time Time("POSIXct" or"POSIXt") for which the PPA needs to be calculated.
#' Use time = c(time1,time2) to calculate PPA for a time range. Default is NULL: calculate PPA for entire STP_track
#' @param points The points used for the PPA calculation given as a vector of integers. Default is NULL: calculate PPA for entire STP_track
#' @param x_density Paramter used for calculating the PPA of entire STPs.
#' The amount of x coordinates for which the corresponding y coordinate(s) will be calculated.
#' Only relevant if the PPA for at least 1 complete STP needs to be calculated.
#' @param time_interval The time interval in minutes used for calculating the PPA.
#' Only used for calculating the PPA for a specfic moment in time and
#' if only a part of the PPA of a STP needs to be calculated.
#' Default is every minute
#' @param quadsegs Passed to buffer. Number of line segments to use to approximate a quarter circle.
#' Only used where paramter time_interval is relavant.
#' @param point_uncertainty uncertainty of space-time points. Default = 0.
#' @return The Potential Path Area as SpatialPolygons.
#' If time is equal to time of a space-time point and point_uncetainty = 0, method returns NA because there is no PPA
#' @author Mark ten Vregelaar
#' @importMethodsFrom raster bind
#' @importFrom rgeos gBuffer gIntersection gUnaryUnion
#' @importFrom rgeos gBuffer gIntersection gUnaryUnion
#' @import plyr rgdal
#' @export
#' @examples
#'library(spacetime)
#'library(sp)
#'#--------------------------create a STP_Track--------------------------
#'# set time
#'t1 <- strptime("01/01/2017 00:00:00", "%m/%d/%Y %H:%M:%S")
#'t2 <- t1+5*60*60
#'time<-seq(t1,t2,30*60)
#'# set coordinates
#'x=c(seq(0,25,5),seq(27.5,37.5,2.5))
#'y=sample(-2:2, 11,replace = TRUE)
#'
#'n = length(x)
#'crs_NL = CRS("+init=epsg:28992")
#'
#'# create class STIDF
#'stidf1 = STIDF(SpatialPoints(cbind(x,y),crs_NL), time, data.frame(co2 = rnorm(n),O2=rnorm(n)))
#'
#'# Track-class {trajectories}
#'my_track1<-Track(stidf1)
#'
#'# set maximum speed
#'v1<-getVmaxtrack(my_track1)+0.00015
#'
#'# STP_track class
#'STP_track1<-STP_Track(my_track1,v1)
#'#------------------------------example 1------------------------------
#'## PPA entire track
#'#calculate PPA
#'PPA<-calculate_PPA(STP_track1)
#'
#'# plot results
#'plot(STP_track1,type='b')
#'plot(PPA,add=TRUE)
#'#------------------------------example 2------------------------------
#'## PPA only using every second point
#'# calculate PPA
#'PPA<-calculate_PPA(STP_track1,points = seq(1,11,2))
#'
#'# plot results
#'plot(STP_track1,type='b')
#'plot(PPA,add=TRUE)
#'#------------------------------example 3------------------------------
#'## PPA of a specfic moment in time
#'# calculate PPA
#'time <- strptime("01/01/2017 01:15:00", "%m/%d/%Y %H:%M:%S")
#'PPA<-calculate_PPA(STP_track1,time = time)
#'
#'# plot results
#'plot(STP_track1,type='b')
#'plot(PPA,add=TRUE)
#'#------------------------------example 4------------------------------
#'## PPA for a time range
#'# calculate PPA
#'timerange1 <- c(t1,strptime("01/01/2017 02:15:00", "%m/%d/%Y %H:%M:%S"))
#'PPA<-calculate_PPA(STP_track1,time = timerange1)
#'
#'# plot results
#'plot(STP_track1,type='b')
#'plot(PPA,add=TRUE)
calculate_PPA <-
  function(STP_track, time = NULL, points = NULL, x_density = 250,
           time_interval = 1, quadsegs = 12, point_uncertainty = 0) {
    if (!is.null(points)) {
      STP_track <- STP_track[points, '']
    }
    if (length(time) == 1) {
      result <-  calc_PPA(STP_track, time[1], qs = quadsegs,point_uncertainty = point_uncertainty)
      if (is.null(result)) {
        stop('Could not calculate PPA. Maximum speed might be to low.')}
      return(result)
      }
     else if (length(time) == 2) {
       if (time[2]<=time[1]){
         stop("error in time. time[2] is smaller or eaual to time[1]")
       }
       else if(time_interval > difftime(time[2],time[1],units = 'mins')){
         stop("error in time. time interval is bigger than time difference between space-time points")

       }else{
         result <- calcPPA_STP_Tinterval(STP_track, time, x_density, time_interval, quadsegs)
       }
       }
         else{
      result <- calcPPA_STP_Track(STP_track, x_density)
    }
    if (!isS4(result)) {
      warning(
        'Could not calculate PPA. Maximum speed might be to low or time is equal to time of a
        space-time point in which case the PPA is a point. Returning NA'
      )
    } else if (point_uncertainty > 0) {
      result <- gBuffer(result, width = point_uncertainty)
    }

    return(result)

  }

# @title calc_PPA
# @description This function calcualtes the Potential Path Area(PPA) for a given moment in time.
# @param STP A STP_Track
# @param time Time("POSIXct" or"POSIXt") for which the PPA needs to be calculated.
# @param points The points used for the PPA calculation given as a vector of integers.
# Default are the space-time points are directlty before and after space-time point
# @return The Potential Path Area as SpatialPolygons for the provided time
# @author Mark ten Vregelaar
# @importFrom rgeos gBuffer gIntersection gUnaryUnion
#
calc_PPA <- function(STP,t,points=NULL,qs=12,point_uncertainty=0){


  # If there are no points provided, find the two points before and after t.
  if (is.null(points)){
    i=1
    search=T
    while (i+1<=length(STP) & search){
      if (t>=STP@endTime[i] & t<=STP@endTime[i+1]){
        p1=i
        p2=i+1
        search = F
      }
      i=i+1

    }
  }
  else{
    p1 = points[1]
    p2 = points[2]
  }

  # get maximum speed
  v = STP@connections$vmax[p1]

  # get time difference of the time between the two points and t
  t1 = abs(difftime(STP@endTime[p1],t,units = 'secs'))
  t2 = abs(difftime(STP@endTime[p2],t,units = 'secs'))

  #calculate the maximum travel dictance strating form the two points in m
  # now assuming meters as unit for the projecition<-----------------------------------------
  s1 = v*as.numeric(t1)
  s2 = v*as.numeric(t2)

  # get coords of start and end point
  startpoint <- STP@sp[p1,]
  endpoint <- STP@sp[p2,]
  # calculate and return the PPA, the intersection of the two buffers
  if(t==STP@endTime[p1] & point_uncertainty>0){
    return(gBuffer(startpoint,width = point_uncertainty))
  }
  if(t==STP@endTime[p2] & point_uncertainty>0){
    return(gBuffer(endpoint,width = point_uncertainty))
  }
  PPA <- tryCatch(
    {
      if (t1>0 & t2 >0){
        # calculate area that can be reached from each point
        buffer1<-gBuffer(startpoint,width=s1,quadsegs=qs)#default 6, first test 7
        buffer2<-gBuffer(endpoint,width=s2,quadsegs=qs)
        PPA<-gIntersection(buffer1,buffer2)}
      else {
        NA}

    },
    error=function(cond) {
      message("No intersection")
      message("Here's the original error message:")
      message(cond)
      # Choose a return value in case of error
      return(NA)
    },
    warning=function(cond) {
      message("functions caused a warning:")
      message("Here's the original warning message:")
      message(cond)
      # return PPA
      return(PPA)
    })
  return(PPA)
}

# @title calcPPA_STP
# @description This function calcualtes the Potential Path Area(PPA) of Space-Time Prism(STP).
# A STP is a set of space-time points with a maximum speed
# @param STP A STP_Track
# @param x_density The amount of x coordinates for which the corresponding y coordinate(s) will be calculated.
# @return The Potential Path Area of the STP as SpatialPolygons
# @author Mark ten Vregelaar
#
calcPPA_STP <- function(STP,x_density=250){
  # get coordiantes of STP
  x1<-STP@sp@coords[1,1]
  y1<-STP@sp@coords[1,2]
  x2<-STP@sp@coords[2,1]
  y2<-STP@sp@coords[2,2]

  # Downsize x and y. For accurate calculation of the PPA the values have to be small
  xmin<-min(x1,x2)
  ymin<-min(y1,y2)
  x1<-x1-xmin
  y1<-y1-ymin
  x2<-x2-xmin
  y2<-y2-ymin


  ## calculate the x coordinates for which the y coordinates need to be calculated
  # total time budget
  t=as.numeric(difftime(STP@endTime[2],STP@endTime[1],units = 'secs'))
  v<-STP@connections$vmax
  # total distance that can be covered
  s=v*t
  # distance that can still be covered after reaching the other point
  dist <- s-STP@connections$distance

  # caclulate values for the x coordinate usingthe remaining distance
  if(x1<x2){
    xrange<-c(x1-dist,x2+dist)
  }else{
    xrange<-c(x2-dist,x1+dist)
  }

  xpoints<-seq(xrange[1],xrange[2],length=x_density)
  # calculate y coordinates. y_segment1 and y_segment2 are both half a ellipse or circle
  suppressWarnings(yCoords<-sapply(xpoints,simplify = "array", function(x0) {
    y_segment1 <- ((s*sqrt(s^4+((-2*y2^2) + 4*y1*y2 - 2*x2^2 + 4*x0*x2 - 2*y1^2 - 2*x1^2 + 4*x0*x1 - 4*x0^2)
                   *s^2 + y2^4 - 4*y1*y2^3 + (2*x2^2 - 4*x0*x2 + 6*y1^2 + 2*x1^2 - 4*x0*x1 + 4*x0^2)
                   *y2^2 + ((-4*y1*x2^2) + 8*x0*y1*x2 - 4*y1^3+((-4*x1^2) + 8*x0*x1 - 8*x0^2)*y1)*y2 +
                     x2^4 - 4*x0*x2^3 + (2*y1^2 - 2*x1^2 + 4*x0*x1 + 4*x0^2)*x2^2 + ((-4*x0*y1^2) + 4*x0*x1^2 - 8*x0^2*x1)
                   *x2 + y1^4 + (2*x1^2 - 4*x0*x1 + 4*x0^2)*y1^2 + x1^4 - 4*x0*x1^3 + 4*x0^2*x1^2)
            + (y2+y1) *s^2 - y2^3 + y1*y2^2 + ((-x2^2) +2*x0*x2 + y1^2 + x1^2 - 2*x0*x1)*y2 + y1*x2^2 - 2*x0*y1*x2 -y1^3 + (2*x0*x1 - x1^2)*y1)/
             (2*s^2 - 2*y2^2 + 4*y1*y2 - 2*y1^2))

    y_segment2 <- ((-(s*sqrt(s^4+((-2*y2^2) + 4*y1*y2 - 2*x2^2 + 4*x0*x2 - 2*y1^2 - 2*x1^2 + 4*x0*x1 - 4*x0^2)
                     *s^2 + y2^4 - 4*y1*y2^3 + (2*x2^2 - 4*x0*x2 + 6*y1^2 + 2*x1^2 - 4*x0*x1 + 4*x0^2)
                     *y2^2 + ((-4*y1*x2^2) + 8*x0*y1*x2 - 4*y1^3+((-4*x1^2) + 8*x0*x1 - 8*x0^2)*y1)*y2 +
                       x2^4 - 4*x0*x2^3 + (2*y1^2 - 2*x1^2 + 4*x0*x1 + 4*x0^2)*x2^2 + ((-4*x0*y1^2) + 4*x0*x1^2 - 8*x0^2*x1)
                     *x2 + y1^4 + (2*x1^2 - 4*x0*x1 + 4*x0^2)*y1^2 + x1^4 - 4*x0*x1^3 + 4*x0^2*x1^2)
              +((-y2)-y1)*s^2 + y2^3 - y1*y2^2 + (x2^2 - 2*x0*x2 - y1^2 - x1^2 + 2*x0*x1)*y2 - y1*x2^2 + 2*x0*y1*x2 + y1^3 + (x1^2 - 2*x0*x1)*y1))/
             (2*s^2 - 2*y2^2 + 4*y1*y2 - 2*y1^2))

    c(y_segment1,y_segment2)
  }))
  # put results in data frame and remove NANs
  coords<-data.frame(x=rep(xpoints,2),y=c(yCoords[1,],yCoords[2,]))
  coords<-coords[complete.cases(coords),]
  # take convex hull
  PPA<-coords[chull(coords),]

  # relocate PPA to orginal position
  PPA$x<-PPA$x+xmin
  PPA$y<-PPA$y+ymin
  # covert PPA to SpatialPolygons class and return result
  Poly <- Polygons(list(Polygon(PPA)),'1')
  PPApoly <- SpatialPolygons(list(Poly),proj4string = STP@sp@proj4string)
  return(PPApoly)
}
# @title calcPPA_STP_Track
# @description This function calcualtes the Potential Path Area(PPA) of an entire STP_track
# @param STP_track The STP_Track for which the PPA will be calculated
# @param x_density The amount of x coordinates for which the corresponding y coordinate(s) will be calculated.
# @return The Potential Path Area of the STP as SpatialPolygons
# @author Mark ten Vregelaar
# @importMethodsFrom raster bind
calcPPA_STP_Track <- function(STP_track,x_density=250){
  if (length(STP_track)==2){

    return(calcPPA_STP(STP_track,x_density))
  }
  ## calculate PPA of entire track
  # Calculate PPA for every STP of the trajectory

  PPAlist<-lapply((1:(length(STP_track)-1)), FUN=function(p1) {

    STP<-STP_track[p1:(p1+1),'']
    # calculate PPA of STP
    calcPPA_STP(STP,x_density)
  })

  #convertlist of PPA to one SpatialPolygons
  PPA_polygons<-do.call(bind,PPAlist)
  PPA<-gUnaryUnion(PPA_polygons)

  #return results
  return(PPA)#<---------------------------------------------------------------------------------------------------------

}

# @title calcPPA_STP_Tinterval
# @description Function for calculating PPA for specific time range
# @param STP_track The STP_track for which the PPA needs to be calculated
# @param time_range time range ("POSIXct" or"POSIXt") for which the PPA needs to be calculated.
# @param x_density Paramter used for calculating the PPA of entire STPs.
# The amount of x coordinates for which the corresponding y coordinate(s) will be calculated.
# Only relevant if the PPA for at least 1 complete STP needs to be calculated
# @param time_interval(minutes) The time interval used for calculating the PPA.
# Only used for calculating the PPA for a specfic moment in time and
# if only a part of the PPA of a STP needs to be calculated.
# Default is every minute
# @param quadsegs Passed to buffer. Number of line segments to use to approximate a quarter circle.
# Only used where paramter time_interval is relavant
# @return The Potential Path Area for the time range as SpatialPolygons
# @author Mark ten Vregelaar
#'
calcPPA_STP_Tinterval <- function(STP_track,time_range,x_density,time_interval, quadsegs){# can be done smarter

  # define frequently used vaiables
  crs<- STP_track@sp@proj4string
  t1<-time_range[1]
  t2<-time_range[2]


  # number of space-time points within time range
  points<-STP_track@endTime>=t1 & STP_track@endTime<=t2
  # if less than two calculate PPAs
  if (sum(points==T)<2){
    # times for which the PPA needs to be calculated
    PPAtimes<-seq(t1+(time_interval*60),t2,(time_interval*60))

      # calculate PPA for times in PPAtimes
      PPAs<-lapply(PPAtimes, function(x) {

        if (!(x %in% STP_track@endTime)){
          calc_PPA(STP_track, x,qs = quadsegs)@polygons[[1]]@Polygons[[1]]@coords
        }
        else{
          NULL# <---- now NULL, alternatives may be faster
        }})
    # test if a space-time point is within the time range
    TF_point <- lapply(STP_track@endTime, FUN=function(x) in_time_range(x,time_range))
    if (T %in% TF_point){
      ## CASE: only one space-time point in time range
      # get times before and after original space-time point
      subset1<-which(PPAtimes<STP_track@endTime[which(TF_point==T)])
      subset2<-which(PPAtimes>STP_track@endTime[which(TF_point==T)])
      # select PPAs for both before and after space-time point
      PPAs1 <-PPAs[subset1]
      PPAs2 <- PPAs[subset2]
      #calculate PPA by calculating convexhull of PPAS
      PPApoly1<-chull_poly(PPAs1,crs)
      PPApoly2<-chull_poly(PPAs2,crs)
      # comind the two PPAs
      PPA_polygons<-do.call(bind,c(PPApoly1,PPApoly2))
      PPA<-gUnaryUnion(PPA_polygons)
      # return final PPA
      return(PPA)
    }

    else{
      ## CASE: no point in between time range
      #calculate PPA by calculating convexhull of PPAS
      PPApoly <- chull_poly(PPAs,crs)
      # return final PPA
      return(PPApoly)
    }}
  else{
    # CASE: at least one complete STP in time-range
    # calculate PPA the complete part of STP_track
    comp_STP_track<-STP_track[1:length(STP_track),paste(time_range[1],time_range[2],sep="::")]
    PPA_track <- calcPPA_STP_Track(comp_STP_track,x_density)
    # initialize PPAs to NULL
    PPApoly1<-NULL
    PPApoly2<-NULL
    ## calc PPA for time before complete part of STP_track

    #lapplay for calculating PPAs for different moments in time
    if ((t1+(time_interval*60))<(comp_STP_track@endTime[1])){
      PPAS1<-lapply(seq(t1+(time_interval*60),comp_STP_track@endTime[1]-(time_interval*60),(time_interval*60)), function(x) {
        calc_PPA(STP_track, x, qs = quadsegs)@polygons[[1]]@Polygons[[1]]@coords
      }
      )
      # convexhull of points
      PPApoly1<-chull_poly(PPAS1,crs)
    }
    ## calc PPA for time after complete part of STP_track

    #lapplay for calculating PPAs for different moments in time
    if ((tail(comp_STP_track@endTime,n=1)+(time_interval*60))<t2){
      PPAS2<-lapply(seq(tail(comp_STP_track@endTime,n=1)+(time_interval*60),t2,(time_interval*60)), function(x)
        calc_PPA(STP_track, x, qs = quadsegs)@polygons[[1]]@Polygons[[1]]@coords)
      # convexhull of points
      PPApoly2<-chull_poly(PPAS2,crs)
    }
    if(is.null(PPApoly1) & is.null(PPApoly2)){
      ## CASE: a complete track
      # no parts before and after complete track
      return(PPA_track)
    }
    else{
      # combine before and after part with complete track
      PPA_polygons<-do.call(bind,c(PPApoly1,PPApoly2,PPA_track))
      PPA<-gUnaryUnion(PPA_polygons)
      return(PPA)
    }
  }
}

# @title in_time_range
# @description Function that checks if the time is within the time_range
# @param time time (POSIXct): time
# @param time_range time_range(two POSIXct): the time range
#
# @return True or False(logical): True if time within time_Range
in_time_range <- function(time, time_range){

  stopifnot(length(time_range) == 2L)

  time>time_range[1] & time<time_range[2]

}

chull_poly <- function(PPAcoords,crs){
  # Function that calculates PPA based on convex hull of PPAcoords
  #
  #   Arg:
  #          PPAcoords(matrix): x and y coordiantes
  #          crs(CRS): crs of coordinates
  #
  #    Return:
  #         PPA(SpatialPolygons): Potential Path Area
  coords<-plyr::rbind.fill.matrix(PPAcoords)
  #coords<-do.call(rbind,PPAS)#<------ other method maybe be faster<________________________
  # convexhull of coords
  coords<-coords[chull(coords),]
  # covert PPA to SpatialPolygons class and return result
  Poly <- Polygons(list(Polygon(coords)),'1')
  PPApoly <- SpatialPolygons(list(Poly),proj4string = crs)
  return(PPApoly)
}



