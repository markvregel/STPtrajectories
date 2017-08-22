 # potential stay

#' @title potential_stay
#' @description Function that calculates the time intervals at which the individual could have been at spgeom,
#'  a spatial location as a point, line or polygon. The total potential stay time might not be equal to the sum of the intervals.
#'  This applies to STP_tracks with multiple intersections for one STP and if the track has time uncertainty.
#'  If the track has time unceratinty the method calculates the max potential stay for each STP,
#'  which may result in a time overlap in the returned potential stay time intervals.
#'  If space-time prism has an activity_time, it is assumed that individual is at the location for the activity.
#'  The potential_stay time interval is thus not affected by the activity_time if the sp_geom can be reached.
#' @param STP_track A STP_Track
#' @param spgeom A Spatialpoints,Spatiallines or Spatialpolygons object
#' @param x_density Paramter used for calculating PPA. Default is 250 coordinates
#' The amount of x coordinates for which the corresponding y coordinate(s) will be calculated.
#' @return A named list with  the potential stay time intervals for each STP that intersects with the spatial geometry
#' @export
#' @importFrom rgeos gIntersects gDistance
#' @examples
#'library(rgeos)
#'library(sp)
#'library(spacetime)
#'library(rgl)
#'#--------------------------create a STP_Track--------------------------
#'# set time
#'t1 <- as.POSIXct(strptime("01/01/2017 12:00:00", "%m/%d/%Y %H:%M:%S"))
#'t2 <- as.POSIXct(strptime("01/01/2017 14:00:00", "%m/%d/%Y %H:%M:%S"))
#'n<-5
#'time<-seq(t1,t2,length.out = n)
#'# set coordinates
#'x=c(1,3,4,7,6)
#'y=c(2,5,8,9,11)
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
#'v1<-getVmaxtrack(my_track1)*1.5
#'
#'# STP_track class
#'STP_track1<-STP_Track(my_track1,v1)
#'#--------------------------create a spatialpolygon--------------------------
#'point_xy<-c(5,9)
#'lake <- readWKT("POLYGON((1 8,1.3 8.5,1.7 8.7,2 9,2.5 9.2,3 9,3.5 8.5,3.2 7.7,2.5 7.5,2 7.5,1 8))")
#'lake@proj4string<-crs_NL
#'
#'# calculate the potential stay time for point
#'intervals <- potential_stay(STP_track = STP_track1,spgeom = lake)
#'
#'# calculate the time the individual could have been at the lake
#'lake_time <- sum(sapply(intervals, function(STP){difftime(STP[2],STP[1],units = 'mins')}))
#'print(paste('Total time individual could have been at the lake is ',round(lake_time,2),'minutes'))
#'
#'# visulise in 2D
#'plot(STP_track1,type='b')
#'PPA_track<-PPA(STP_track1)
#'plot(PPA_track,add=T)
#'plot(lake,add=T,border= 'blue',lwd=2)
#'
#'# visulise in 3D
#'open3d()
#'zfac<-STP_plot(STP_track1,time_interval = 0.6)
#'
#'axes_STP_plot(c(STP_track1@endTime[1],STP_track1@endTime[n]),z_factor = zfac,n_ticks_z = 5,n_ticks_xy = 4)
#'
#'x<-lake@polygons[[1]]@Polygons[[1]]@coords[,1]
#'y<-lake@polygons[[1]]@Polygons[[1]]@coords[,2]
#'z<-difftime(STP_track1@endTime[n],STP_track1@endTime[1],units = 'mins')*zfac
#'
#'shade3d(translate3d(
#'  extrude3d(x,y,thickness = z),0,0,0),col='blue',add=TRUE)
potential_stay <- function(STP_track, spgeom, x_density = 250) {
  # calculate PPA
  PPA_track <- PPA(STP_track,x_density = x_density)
  # initialise intervals to store time intervals of potential stay
  intervals <- c()
  # check if there is a potential stay time
  if (gIntersects(PPA_track, spgeom)) {
    # time uncertainty in seconds
    tu <- STP_track@rough_sets$time_uncertainty*60


    # check which STPs intersect with the spgeom
    for (i in 1:length(STP_track@connections[, 1])) {
      STP <- STP_track[i:(i + 1), '']
      PPA_STP <- PPA(STP,x_density = x_density)

      if (gIntersects(PPA_STP, spgeom)) {
        # if STP intersect with spgeom get intersection
        intersection <- gIntersection(PPA_STP,spgeom)

        if(!is(intersection,"SpatialPoints")){
          intersection <- disaggregate(intersection)

        if (length(intersection)>1){
          # if intersection results in more than 1 polygon or line then there are multiple time intervals
          warning(paste0("More than one interval for STP " , i,
                        ". Total potential stay time at input location is not equal to
                        the sum of time difference between the intervals"))

        }}

          for (j in 1:length(intersection)){
        # calcualte minimal distance to reach intersection for both points
        dist1 <- gDistance(STP@sp[1,],intersection[j])-STP@rough_sets$location_uncertainty
        dist2 <- gDistance(STP@sp[2,],intersection[j])-STP@rough_sets$location_uncertainty

        if( dist1 < 0)
          dist1<- 0
        if (dist2 < 0)
          dist2<- 0

        # calculate the time the individual can reach closest point of spgeom from  space-time point
        time1 <- STP@endTime[1]-tu+dist1/STP@connections$vmax
        time2 <- STP@endTime[2]+tu-dist2/STP@connections$vmax

        # list intervals and add stp number as name
        STP_int <- list(c(time1,time2))
        names(STP_int)<-paste("STP",i)
        intervals <- c(intervals, STP_int)
      }}

    }
  }
  # if there is not potential stay time, give a warning.
  else{
    warning('No potential stay time. Indivivdual cannot reach spgeom.')
  }
  # return list sorted on time
  return(intervals[order(sapply(intervals,'[[',1))])
}


