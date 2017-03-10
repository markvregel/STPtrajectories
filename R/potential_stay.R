# potential stay

#' @title potential_stay
#' @description
#' @param potential_stay
#' @param spgeom point,line or polygon object
#' @param time_interval
#'
#' @return
#' @export
#'
#' @examples
potential_stay <- function(STP_track, spgeom, time_interval) {
  time_interval<-time_interval * 60
  PPA_track <- calculate_PPA(STP_track)
  intervals <- c()

  if (gIntersects(PPA_track, spgeom)) {
    for (i in 1:length(STP_track@connections[, 1])) {
      STP <- STP_track[i:(i + 1), '']
      PPA <- calculate_PPA(STP)

      if (gIntersects(PPA, spgeom)) {
        times <-seq(STP_track@endTime[i] + time_interval*0.01,
            STP_track@endTime[i + 1] - time_interval*0.01,
            time_interval  )
        good_times <- lapply( 1:length(times),FUN = function(i) {
            PPA <- calculate_PPA(STP, times[i])

            if (gIntersects(PPA, spgeom)) {
              times[i]
            }
          }
        )
        good_times[sapply(good_times, is.null)] <- NULL

        minmax_time <- list(c(good_times[[1]], good_times[length(good_times)]))
        names(minmax_time)<-paste("STP",i)
        intervals <- c(intervals, minmax_time)
      }

    }
  }
  return(intervals)
}


# # points
# STP_track<-STP_moving
# point<-c(1499000,6812200)
# point<-SpatialPoints(matrix(point,nrow = 1),STP_track@sp@proj4string)
# times<-potential_stay(STP_track,spgeom = point,time_interval = 1)
# point2<-STP_track@sp[8]
# times<-potential_stay(STP_track,spgeom = point2,time_interval = 1/60)
#
#
# STP_plot(STP_track,0.5,zfactor = 100,col='blue')
# axes_STP_plot(c(STP_track@endTime[1],STP_track@endTime[length(STP_track)]),z_factor = 100,n_ticks_z = 10)
# point_buffer<-gBuffer(point2,width = 50)
# x<-point_buffer@polygons[[1]]@Polygons[[1]]@coords[,1]
# y<-point_buffer@polygons[[1]]@Polygons[[1]]@coords[,2]
# shade3d(translate3d(
#   extrude3d(x,y,thickness = 25000),0,0,0),col='red',add=TRUE)
#
# # polygons
#
# lake<-gBuffer(STP_track@sp[4],width = 200)
# times<-potential_stay(STP_track,spgeom = lake,time_interval = 60)
#
#
# open3d()
# STP_plot(STP_track,0.5,zfactor = 100,col='red')
# axes_STP_plot(c(STP_track@endTime[1],STP_track@endTime[length(STP_track)]),z_factor = 100,n_ticks_z = 10)
# x<-lake@polygons[[1]]@Polygons[[1]]@coords[,1]
# y<-lake@polygons[[1]]@Polygons[[1]]@coords[,2]
# shade3d(translate3d(
#   extrude3d(x,y,thickness = 25000),0,0,0),col='blue',add=TRUE)
#
# # experimenting with optimizers
#
#
# ## "wrong" solution with unlucky interval and piecewise constant f():
# f  <- function(x) ifelse(x > -1, ifelse(x < 4, exp(-1/abs(x - 1)), 10), 10)
# fp <- function(x) { print(x); f(x) }
#
# plot(f, -2,5, ylim = 0:1, col = 2)
# optimize(fp, c(-4, 20))   # doesn't see the minimum
# optimize(fp, c(-7, 20))   # ok
#
# f <- function(x){
#   print(as.POSIXct(x,origin = "1969-12-31 23:48:43",tz = "CET"))
#   x=as.POSIXct(x,origin = "1970-01-01",tz = "CET")
#   PPA<-calculate_PPA(STP,x)
#   if(gIntersects(PPA,point))
#   {0
#   }else{
#       1}}
#
# optimize(f,c(STP@endTime[1],STP@endTime[2]),maximum = F)
# optim(f,c(STP@endTime[1],STP@endTime[2]))
# x=STP@endTime[1]
# strptime(x)
#
#
#
# fr <- function(x) {   ## Rosenbrock Banana function
#   x1 <- x[1]
#   x2 <- x[2]
#   100 * (x2 - x1 * x1)^2 + (1 - x1)^2
#   print(100 * (x2 - x1 * x1)^2 + (1 - x1)^2)
# }
# grr <- function(x) { ## Gradient of 'fr'
#   x1 <- x[1]
#   x2 <- x[2]
#   c(-400 * x1 * (x2 - x1 * x1) - 2 * (1 - x1),
#     200 *      (x2 - x1 * x1))
# }
# optim(c(-1.2,1), fr)
