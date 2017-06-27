#visulize STPs
# fix Cannot triangulate polygon
#' @title STP_plot
#' @description This function visualizes STPs in 3D.
#' @param STP_track A STP_Track.
#' @param time_interval the time interval in minutes.Determines the amount of PPAs that are plotted.
#' @param zfactor realtive size of z axis compared to x and y axis. default=NULL: ratio axes depns on input data
#' @param color of STP(s).
#' @param st start time as "POSIXct" or"POSIXt". For plotting multiple STP_tracks use 1 starttime.
#' @param alpha transparency of STPs. between 0 and 1
#' @param cut_prisms In case of time uncerainty: wheter to cut the middle and
#' normally overlapping prisms at the original time of the space-time points.
#' Only relevant if alpha<1. default is TRUE
#' @param quadsegs Passed to PPA Number of line segments to use to approximate a quarter circle.
#' @importFrom  rgl translate3d extrude3d shade3d rgl.triangles
#' @importFrom geometry convhulln
#' @author Mark ten Vregelaar
#' @return If no zfactor is provided, method returns calculated zfactor.
#' Otherwise mehtod returns NULL
#' @export
#' @examples
#'library(spacetime)
#'library(sp)
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
#'
#'## 3D STP plot of STP_tracks
#'z_fac<-0.2 # relative size of z scale/aspect ratio to spatial scale
#'# plot STPS first STP_track
#'STP_plot(STP_track1,time_interval = 1,z_fac)
#'# plot STPS second STP_track
#'STP_plot(STP_track2,time_interval = 1,z_fac,'blue',st = STP_track1@endTime[1])
#'# provide st for correct starting location first STP
#'
#'# calculate first and last moment in time
#'min_max_Time<-c(STP_track1@endTime[1],STP_track2@endTime[length(STP_track2)])
#'# add axes
#'axes_STP_plot(min_max_Time,z_factor = z_fac)
#'
#'# add title and change background colour
#'library(rgl)
#'title3d(main = '2 randomly generated STP tracks')
#'bg3d('lightblue')
STP_plot<-function(STP_track,time_interval=0.5,zfactor=NULL,col='red',
                   st=NULL,alpha=1,cut_prisms=TRUE,quadsegs=12){
  # get length of STP_track
  n <- length(STP_track)
  # time_interval in seconds
  time_interval<- time_interval*60

  # time uncertainty in secs
  tu <- STP_track@rough_sets$time_uncertainty*60

  # if no zfacor provided calculate one based on longets spatial axis
  if(is.null(zfactor)){
    bbox<-bbox(STP_track)
    spatial_diff<-max(bbox[1,2]-bbox[1,1],bbox[2,2]-bbox[2,1])
    time_diff<-difftime(STP_track@endTime[[n]]+tu,STP_track@endTime[[1]]-tu,units = 'mins')
    zfac<-spatial_diff/as.numeric(time_diff)

    if (zfac==0){
      stop('could not calculate z factor because first and last point have same location.
           retry with value for zfactor paramter')
    }
  }else{
  zfac<-zfactor
}

  # if no st provided. start at first time of STP_track
  if (is.null(st)){
    st <- STP_track@endTime[[1]]-tu
  }else{
    if (length(st)!=1){
      stop('startime is not valid. st should have length 1')

    }
  }


  # loop through STPS
  for(i in 1:(length(STP_track)-1)){

  # create list of times for which PPAS need to be calculated
  times<-seq(STP_track@endTime[i],STP_track@endTime[i+1],(time_interval))
  if(tail(times,1)!=STP_track@endTime[i+1]){
  times<-c(times,STP_track@endTime[i+1])
  }
  at <- STP_track@connections$activity_time[i]
  # if there is activty time calculate normal PPA and clip to smaller PPA with at
  if(at>0){
    STP <- STP_track[i:(i+1),'']
    PPA_STP <-PPA(STP)
    STP_track@connections$activity_time[i]<-0
    suppressWarnings(PPAS<-lapply(times, function(x) {
      PPA1<-PPA(STP_track,points = c(i,i+1), x,quadsegs = quadsegs)
      # only intersect if PPA could be calcualted
      if(is.na(PPA1)){
        return(PPA1)
        }else{
          gIntersection(PPA1,PPA_STP)
        }
    }))
    STP_track@connections$activity_time[i]<-at
  }else{

  # calculate PPAS
  suppressWarnings(PPAS<-lapply(times, function(x) {
    PPA(STP_track, x,points = c(i,i+1),quadsegs = quadsegs)
  }))
    }

  #create times and ppas bottum and top of STP_track in case of time uncertainty
  # and if cut_prisms is FALSE or if stp is first or last
   if(tu>0 & (i==1 | i==n-1 | !cut_prisms)){

    STP1 <- STP_track[i:(i+1),'']
    STP1@rough_sets$time_uncertainty <- 0
    STP1@endTime<-c(STP1@endTime[1] - tu,STP1@endTime[2] + tu)
    #bottum
    if(i==1 | !cut_prisms){
      times1 <- seq(STP_track@endTime[i]-tu,STP_track@endTime[i]-time_interval,time_interval)
      suppressWarnings(PPAS1<-lapply(times1, function(x) {
        PPA(STP1, x)
      }))
      PPAS<-c(PPAS1,PPAS)
      times<-c(times1,times)
    }
    #top
  if(i==(n-1) | !cut_prisms){

    times2 <- seq(STP_track@endTime[i+1]+time_interval,STP_track@endTime[i+1]+tu,time_interval)
    suppressWarnings(PPAS2<-lapply(times2, function(x) {
      PPA(STP1, x)
    }))
    PPAS<-c(PPAS,PPAS2)
    times<-c(times,times2)


  }
  }


  # remove NAs for PPAs that could not be calculated
  NAs<-!is.na(PPAS)

  PPAS<-PPAS[NAs]

  # nummerical list of times for which PPAS were calculated
  tdif<-as.numeric(difftime(times,st,units = 'min'))*zfac
  t<-tdif[NAs]

  # get xyz coordinates for plotting
  STP_coords<-lapply(1:length(PPAS),function(j){
    coords<-PPAS[[j]]@polygons[[1]]@Polygons[[1]]@coords
    x<-coords[,1]
    y<-coords[,2]
    z<-rep(t[j],length(x))
    list(x=x,y=y,z=z)
    })
  STP_coords<-do.call(Map,c(c,STP_coords))

  # add orignal space-time points to STP if there is no point uncerainty.
  # results in overlapping STPs if tu >0 & cut_prisms==FAlSE
  if(STP_track@rough_sets$location_uncertainty==0){
    # add 1st control point
    if(!cut_prisms | tu ==0 | i==1){
    STP_coords$x<-c(STP_track@sp@coords[i,1],STP_coords$x)
    STP_coords$y<-c(STP_track@sp@coords[i,2],STP_coords$y)
    STP_coords$z<-c((as.numeric(difftime((STP_track@endTime[(i)]-tu),st,units = 'min')))*zfac,STP_coords$z)
    }
    # add 2nd control point
    if(!cut_prisms | tu ==0 | i==(n-1)){
    STP_coords$x<-c(STP_coords$x,STP_track@sp@coords[(i+1),1])
    STP_coords$y<-c(STP_coords$y,STP_track@sp@coords[(i+1),2])
    STP_coords$z<-c(STP_coords$z,(as.numeric(difftime((STP_track@endTime[(i+1)]+tu),st,units = 'min')))*zfac)
    }
  }

  # matrix with all coordinates
  stp3d<-do.call(cbind,STP_coords)



    # take convexhull and plot STP
    conv<-t(convhulln(stp3d))
    rgl.triangles(stp3d[conv,1],stp3d[conv,2],stp3d[conv,3],col=col,alpha=alpha)


    }
  # return zfactor is it was not provided
  if(is.null(zfactor)){
return(zfac)}

}

#' @title axes_STP_plot
#' @description This function adds a bbox with axis to a STP_plot
#' @param minmaxT a vector of length 2 with two "POSIXct" or"POSIXt" values.
#' The first and last moment in time of plotted tracks.
#' Make sure first time is equal to the tracks that starts as first.
#' Also take into account time_uncetrainty.
#' @param z_factor the z facfor used in the plot
#' @param n_ticks_xy number of ticks used for the x and y axes
#' @param n_ticks_z number of ticks used for the z axes
#' @param expand for expanding the bbox. passed to rgl.bbox.
#' If not all z/time tick are visisble, increase expand value.
#' @importFrom rgl rgl.bbox axes3d
#' @author Mark ten Vregelaar
#' @export
#'
axes_STP_plot<-function(minmaxT,z_factor,n_ticks_xy=3,n_ticks_z=5,expand=1.1){
  tdif<-as.numeric(difftime(minmaxT[2],minmaxT[1],units = 'mins'))
  tickval<-seq(0,tdif*z_factor,length.out = n_ticks_z)
  timesval<-seq(minmaxT[1],minmaxT[2],length.out = n_ticks_z)


    rgl.bbox(xat = 0,yat = 0,zat = tickval, zlab = timesval,color = c("black", "black"), emission = "#252526",
           specular = "#636363", shininess = 5, alpha = 0.8,expand=expand,ylab = 'y)',xlab = 'x')
  axes3d(c('x--', 'x+-', 'y--', 'y+-'),nticks = n_ticks_xy)

}
