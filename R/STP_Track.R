#' @title STP_Track class
#' @description A class to represent Space-Time Prism(STP) trajecories.
#' These are trajectories with a maximum speed for each segment.
#' The maximum speed is added to the connections slot
#' @param Track object of class  \link{Track}
#' @inheritSection trajectories::Track Slots of class "Track"
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
  if (!all(diff(i)==1)){
    warning('Some intermediate space-time points will be removed. If vmax values differ along the STP_track, the vmax values might change.
The maxmimum speed of new line segments will depend on the first point of the segment')
  }
  track<-Track(as(x, "STIDF")[i, j, ..., drop = drop])

  vmax<-x[['vmax']][x@endTime %in% track@endTime]
  STP_Track(track,head(vmax,-1))

}

setMethod("[", "STP_Track", subs.STP_Track)
