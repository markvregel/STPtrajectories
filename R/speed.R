#' @title getMinimalSpeed
#' @description   This functions calculates the minimal speed required to reach every point
#' The  speed is based on linear movement between two points.
#' @param track the trajectory as \link{STP_Track} or \link{Track}
#' @return Vector of minimal speeds in unit spatial projection unit/s The speed required to reach the next point
#' in the available time
#' @author Mark ten Vregelaar
#' @export
getMinimalSpeed<- function(track){
  # get the number of points of the trajectory
  n<-length(track)

  # calculate the time difference between consecutive points
  tdif <- difftime(track@endTime[2:n], track@endTime[1:(n-1)],
                   units="secs")
  tdif <- as.numeric(tdif)

  # calculate the required speed to reach every point
  minSpeeds <- track@connections$distance/tdif
  return(minSpeeds)
}

#' @title getVmaxtrack
#' @description This functions calculates the maximum speed found in a trajectory.
#' The maximum speed is based on linear movement between measured points.
#' @param track the trajectory as \link{STP_Track} or \link{Track}
#' @return max speed of the moving object
#' @author Mark ten Vregelaar
#' @export
getVmaxtrack <- function(track){
  # calculate the required speed to reach every point
  speeds<- getMinimalSpeed(track)
  # return the maximum speed
  return(max(speeds))
}
