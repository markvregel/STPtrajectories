#' @title alibi_query
#' @description CONTAINS ERRORS. PROBLEMS WITH CASE 3!!!!!
#' This function tests wether there was a possible meeting between two individuals or other moving objects.
#' If the individuals of two trajecories could have met is tested by applying the alibi query to segments that overlap in time.
#' The alibi query is a Boolean query that checks whether two moving individuals,
#' that are given by two samples of space-time points and speed limitations, could have met each other.
#' The query tests if two space-time prisms intersect.
#' Kuijpers et al. (2011) provide the analytical solution for the alibi query that is used by this function
#' @param STP_track1 STP_track1
#' @param STP_track2 STP_track2
#' @param stop_if_true logigal:Stop if intersection is found. Default=TRUE
#' @return If TRUE returns vector space-time point of intersecting STPs. IF no intersection is found returns FALSE
#' @author Mark ten Vregelaar
#' @references - 	Kuijpers, B., Grimson, R., & Othman, W. (2011).
#' An analytic solution to the alibi query in the space-time prisms model for moving object data.
#' International Journal of Geographical Information Science, 25(2), 293-322.
#' @importFrom lubridate interval int_overlaps
#' @importFrom rgeos gConvexHull
#' @export
#' @examples
#'library(spacetime)
#'library(sp)
#'
#'## create 2 STP_tracks
#'# time
#'t1 <- strptime("01/01/2017 00:00:00", "%m/%d/%Y %H:%M:%S")
#'t2 <- t1+5*60*60 # 5 hours after t1
#'time1<-seq(t1,t2,30*60)
#'time2<-time1+0.25*60*60
#'# spatial coordinates
#'x1=c(seq(0,25,5),seq(27.5,37.5,2.5))
#'y1=sample(-2:2, 11,replace = TRUE)
#'x2=c(seq(0,25,5),seq(27.5,37.5,2.5))
#'y2=sample(-2:2, 11,replace = TRUE)
#'
#'n = length(x1)
#'crs_NL = CRS("+init=epsg:28992")
#'
#'# create class STIDF
#'stidf1 = STIDF(SpatialPoints(cbind(x1,y1),crs_NL), time1, data.frame(co2 = rnorm(n),O2=rnorm(n)))
#'stidf2 = STIDF(SpatialPoints(cbind(x2,y2),crs_NL), time2, data.frame(co2 = rnorm(n),O2=rnorm(n)))
#'
#'# Track-class {trajectories}
#'my_track1<-Track(stidf1)
#'my_track2<-Track(stidf2)
#'# set maximum speed
#'v1<-getVmaxtrack(my_track1)+0.00015
#'v2<-getVmaxtrack(my_track2)+0.00030
#'# STP_track class
#'STP_track1<-STP_Track(my_track1,v1)
#'STP_track2<-STP_Track(my_track2,v2)
#'
#'## the alibi query
#'alibi_query(STP_track1,STP_track2)
alibi_query<-function(STP_track1,STP_track2,stop_if_true=TRUE,return_PIA=FALSE,time_interval=1){

  # no rough_sets for alibi query
  STP_track1<-zero_rough_sets(STP_track1)
  STP_track2<-zero_rough_sets(STP_track2)

  Switch<-F
# determine 1st and 2nd STP_track to improve processing time
if (length(STP_track1)!=length(STP_track2)){

  if (length(STP_track1)>length(STP_track2)){
    track1<-STP_track1
    track2<-STP_track2
  }

  if (length(STP_track1)<length(STP_track2)){
    track1<-STP_track2
    track2<-STP_track1
    Switch<-T
  }
} else if(difftime(max(STP_track1@endTime),min(STP_track1@endTime))>difftime(max(STP_track2@endTime),min(STP_track2@endTime))){
  track1<-STP_track1
  track2<-STP_track2
}else{
  track1<-STP_track2
  track2<-STP_track1
  Switch<-T
}


  track2_int<-interval(min(track2@endTime),max(track2@endTime))
  Trues<-c()
#apply alibi query to segments in same time interval
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

          if(Switch){
            STPs<-c(j,i)
          }else{
            STPs<-c(i,j)
          }
          #message(c(TRUE, paste(': possibe intersection for STPs/connections ',STPs[1],'and',STPs[2])))
          points<-list('a'=STPs[1]:(STPs[1]+1),
            'b'=STPs[2]:(STPs[2]+1))

          names(points)<-c('STP_track1','STP_track2')
          if (return_PIA){
            # calculate PIA and time-interval meeting
            time_PIA<-calc_PIA(STP1,STP2,time_interval)
            points<-c(points,time_PIA)}
          if(stop_if_true){
            return(points)
          }else{
            Trues <- c(Trues,list(points))
          }
        }
      }
    }
  }

}
if (length(Trues)>0){
  return(Trues)
  }else{
  return(FALSE)
}
}


alibi_STP <- function(STP1,STP2){

  # unpack track values


  ts<-c(STP1@endTime[1],STP1@endTime[2],STP2@endTime[1],STP2@endTime[2])
  tmin<-min(ts)
  ts<-as.numeric(difftime(ts,tmin,units = 'secs'))#<--------------------------
  t1<-ts[1];t2<-ts[2];t3<-ts[3];t4<-ts[4]

  x1 <- STP1@sp@coords[1,1]
  x2 <- STP1@sp@coords[2,1]

  y1 <- STP1@sp@coords[1,2]
  y2 <- STP1@sp@coords[2,2]

  v1 <- STP1@connections$vmax

  x3 <- STP2@sp@coords[1,1]
  x4 <- STP2@sp@coords[2,1]

  y3 <- STP2@sp@coords[1,2]
  y4 <- STP2@sp@coords[2,2]

  v2 <- STP2@connections$vmax

  result<-alibi(t1, x1, y1, t2, x2, y2, v1, t3, x3, y3, t4, x4, y4,v2)

  query_result<-FALSE
  if(is.numeric(result)){
    #message(paste0("Case ",result))
    if(result==3){
      if(gIntersects(PPA(STP1),PPA(STP2))){
        query_result<-TRUE
      }else{
      }}else{
        query_result<-TRUE
      }}




  return(query_result)
}

calc_PIA<-function(STP1,STP2,time_interval){

  # calculate potential PIA
  PPA1<-PPA(STP1)
  PPA2<-PPA(STP2)
  PPIA<-gIntersection(PPA1,PPA2)

  # Calculate time interval for which meeting is possible
  STP1_time<-potential_stay(STP1,PPIA)[[1]]
  STP2_time<- potential_stay(STP2,PPIA)[[1]]

  if(STP1_time[1]>STP2_time[1]){
    st <- STP1_time[1]
  }else{
    st <- STP2_time[1]
  }

  if(STP1_time[2]<STP2_time[2]){
    et <- STP1_time[2]
  }else{
    et <- STP2_time[2]
  }
  # remove sure time control points not in time because cannot calculate PPA for the control points
  if(et == STP1@endTime[2] | et == STP2@endTime[2] ){
    et<-et-time_interval
  }

  if(st == STP1@endTime[1] | st == STP2@endTime[1] ){
    st<-st+time_interval
  }
  if (st>et){
    return(list("False positive by alibi query. No meeting possible"))
  }
  times<-seq(st,et,time_interval)# times for which will be tested if intersection is possible
  PIAs<-c()
  t1_unknown=T
  for (i in 1:length(times)){
    t <- times[i]


    PPA_STP1<-PPA(STP1,t)
    PPA_STP2<-PPA(STP2,t)


    if(gIntersects(PPA_STP1,PPA_STP2)){
      if(t1_unknown){
        t1i<-i
        t1_unknown=FALSE}

      PIA<-gIntersection(PPA_STP1,PPA_STP2)
      PIAs<-c(PIAs,PIA)

    }
  }

  if (t1_unknown){
    return(list("False positive by alibi query. No meeting possible"))
  }else{
    t1<-times[t1i]
    t2<- times[t1i+length(PIAs)-1]
  PIA_polygons<-do.call(bind,PIAs)
  PIA<-gConvexHull(PIA_polygons)
  # if case 1: control point in PIA, adjust t1 and t2
  if(gIntersects(STP1@sp[1,],PIA) | gIntersects(STP2@sp[1,],PIA)){
    t1<-t1-time_interval
    }

  if(gIntersects(STP1@sp[2,],PIA) | gIntersects(STP2@sp[2,],PIA)){
    t2<-t2+time_interval}

  return(list(meeting_time=c(t1,t2),PIA=PIA))
}}

zero_rough_sets<-function(STP_track1){
  STP_track1@rough_sets$location_uncertainty<-0
  STP_track1@rough_sets$time_uncertainty<-0
  return(STP_track1)
}

