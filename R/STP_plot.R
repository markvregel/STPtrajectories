#visulize STPs
# fix Cannot triangulate polygon
#' @title STP_plot
#' @description This function visualizes STPs in 3D
#' @param STP_track A STP_Track
#' @param time_interval the time interval in minutes.Determines the amount of PPAs that are plotted
#' @param zfactor realtive size of z axis compared to x and y axis
#' @param color of STP(s)
#' @param st start time as "POSIXct" or"POSIXt". For plotting multiple STP_tracks use 1 starttime
#' @importFrom  rgl translate3d extrude3d shade3d
#' @author Mark ten Vregelaar
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
STP_plot<-function(STP_track,time_interval,zfactor=1,col='red',st=NULL){
  #   Return:
  #         No return
  # if no st provided. start at first time of STP_track
  if (is.null(st)){
    st <- STP_track@endTime[[1]]
  }

  # get length of STP_track
  n <- length(STP_track)

  #create list of times for which PPAS need to be calculated
  times<-seq(STP_track@endTime[[1]],STP_track@endTime[[n]]-time_interval*0.01,(time_interval*60))

  #calculate PPAS
  PPAS<-lapply(times, function(x) {
    calc_PPA(STP_track, x)
  })
  # remove NAs for PPAs that could not be calculated
  NAs<-!is.na(PPAS)
  PPAS<-PPAS[NAs]


  # nummerical list of times for which PPAS were calculated
  t<-as.numeric(difftime(times,st,units = 'min'))*zfactor
  t<-t[NAs]

  # plot PPAS in loop lapplay---------------------------------------------------------------
  for (i in 1:length(PPAS)){
    x<-PPAS[[i]]@polygons[[1]]@Polygons[[1]]@coords[,1]
    y<-PPAS[[i]]@polygons[[1]]@Polygons[[1]]@coords[,2]
    tryCatch({
      shade3d(translate3d(extrude3d(x,y,thickness = time_interval*zfactor),0,0,t[i]),col=col,add=TRUE)

    },error=function(cond) {
      message("Whoops could not plot polygon")
      message("Here's the original error message:")
      message(cond)
      # Choose a return value in case of error
      return(NA)})

  }

}

#' @title axes_STP_plot
#' @description This function adds a bbox with axis to a STP_plot
#' @param minmaxT a vector of length 2 with two "POSIXct" or"POSIXt" values
#' @param z_factor the z facfor used in the plot
#' @param n_ticks_xy number of ticks used for the x and y axes
#' @param n_ticks_z number of ticks used for the z axes
#' @importFrom rgl rgl.bbox axes3d
#' @author Mark ten Vregelaar
#' @export
#'
axes_STP_plot<-function(minmaxT,z_factor,n_ticks_xy=3,n_ticks_z=5){
  tdif<-as.numeric(difftime(minmaxT[2],minmaxT[1],units = 'mins'))
  tickval<-seq(0,tdif*z_factor,length.out = n_ticks_z)

  timesval<-seq(minmaxT[1],minmaxT[2],length.out = n_ticks_z)

  rgl.bbox(xat = 0,yat = 0,zat = tickval, zlab = timesval,color = c("black", "black"), emission = "#252526",
           specular = "#636363", shininess = 5, alpha = 0.8,expand=1.02,ylab = 'y)',xlab = 'x')
  axes3d(c('x--', 'x+-', 'y--', 'y+-'),nticks = n_ticks_xy)

}

