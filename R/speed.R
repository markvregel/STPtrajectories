#' @title getVmaxtrack
#' @description This functions calculates the maximum speed found in a trajectory.
#' The maximum speed is based on linear movement between measured points.
#' @param track the trajectory as \link{STP_Track} or \link{Track}
#' @return max speed of the moving object
#' @author Mark ten Vregelaar
#' @export
#' @examples
#'## get a track dataset
#'data(A3)# see trajectories package: importenvirocar
#'car_track<-Track(A3)# recalculate connections data
#'head(car_track@connections)# distance in meters and speed in m/s
#'## get maximum speed
#'speed_ms <- getVmaxtrack(car_track)
#'# maximum speed. It is the speed required to be able to reach every point in time in km/h
#'speed_ms*3.6
getVmaxtrack <- function(track){
  # return the maximum speed
  return(max(track@connections$speed))
}
