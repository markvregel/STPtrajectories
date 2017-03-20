# potential stay

#' @title potential_stay
#' @description
#' @param potential_stay
#' @param spgeom point,line or polygon object
#' @param time_interval
#' @return
#' @export
#' @importFrom rgeos gIntersects
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
        times <-seq(STP@endTime[1] + time_interval*0.01,
                    STP@endTime[2] - time_interval*0.01,
            time_interval)
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


