
alibi_STP <- function(STP1,STP2){


  ts<-c(STP1@endTime[1],STP1@endTime[2],STP2@endTime[1],STP2@endTime[2])
  tmin<-min(ts)
  ts<-abs(as.numeric(difftime(tmin,ts,units = 'hours')))#<--------------------------
  t1<-ts[1];t2<-ts[2];t3<-ts[3];t4<-ts[4]

  x1 <- STP1@sp@coords[1,1]
  x2 <- STP1@sp@coords[2,1]

  y1 <- STP1@sp@coords[1,2]
  y2 <- STP1@sp@coords[2,2]

  v1 <- STP1@connections$vmax*3600

  x3 <- STP2@sp@coords[1,1]
  x4 <- STP2@sp@coords[2,1]

  y3 <- STP2@sp@coords[1,2]
  y4 <- STP2@sp@coords[2,2]

  v2 <- STP2@connections$vmax*3600
  # print(c(t1,t2,t3,t4))
  # print(c(t1, x1, y1, t2, x2, y2, v1))
  # print(c(t3, x3, y3, t4, x4, y4,v2))
  return(alibi(t1, x1, y1, t2, x2, y2, v1, t3, x3, y3, t4, x4, y4,v2))
}


#' @title alibi_query
#' @description CONTAINS ERRORS. CASE 3 goed wrong!!!
#' This function tests wether there was a possible meeting between two individuals or other moving objects.
#' If the individuals of two trajecories could have met is tested by applying the alibi query to segments that overlap in time.
#' The alibi query is a Boolean query that checks whether two moving individuals,
#' that are given by two samples of space-time points and speed limitations, could have met each other.
#' The query tests if two space-time prisms intersect.
#' Kuijpers et al. (2011) provide the analytical solution for the alibi query that is used by this function
#' @param STP_track1 STP_track1
#' @param STP_track2 STP_track2
#' @return True or False for the alibi query
#' @author Mark ten Vregelaar
#' @references - 	Kuijpers, B., Grimson, R., & Othman, W. (2011).
#' An analytic solution to the alibi query in the space-time prisms model for moving object data.
#' International Journal of Geographical Information Science, 25(2), 293-322.
#' @importFrom lubridate interval int_overlaps
#' @export
alibi_query<-function(STP_track1,STP_track2){

# determine 1st and 2nd STP_track to improve processing time
if (length(STP_track1)!=length(STP_track2)){

  if (length(STP_track1)>length(STP_track2)){
    track1<-STP_track1
    track2<-STP_track2
  }

  if (length(STP_track1)<length(STP_track2)){
    track1<-STP_track2
    track2<-STP_track1
  }
} else if(difftime(max(STP_track1@endTime),min(STP_track1@endTime))>difftime(max(STP_track2@endTime),min(STP_track2@endTime))){
  track1<-STP_track1
  track2<-STP_track2
}else{
  track1<-STP_track2
  track2<-STP_track1
}
  track2_int<-interval(min(track2@endTime),max(track2@endTime))
#apply alibi query to segments is same time interval
for (i in 1:(length(track1)-1)){
  stp1_int<-interval(track1@endTime[i],track1@endTime[i+1])
  if(int_overlaps(stp1_int,track2_int)){

    for (j in 1:(length(track2)-1)){
      stp2_int<- interval(track2@endTime[j],track2@endTime[j+1])
      if (int_overlaps(stp1_int,stp2_int)){
        STP1 <- track1[i:(i+1),'']
        STP2 <- track2[j:(j+1),'']
        result<-alibi_STP(STP1,STP2)

        if (result){
          print(c(i,j, "True"))
          return(result)

        }
      }
    }
  }

}

return(FALSE)

}





