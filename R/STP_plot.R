#visulize STPs
# fix Cannot triangulate polygon
#' @title STP_plot
#' @description This function visualizes STPs in 3D.
#' @param STP_track A STP_Track.
#' @param time_interval the time interval in minutes.Determines the amount of PPAs that are plotted.
#' @param zfactor realtive size of z axis compared to x and y axis. default=NULL: ratio axes depns on input data
#' @param color of STP(s).
#' @param st start time as "POSIXct" or"POSIXt". For plotting multiple STP_tracks use 1 starttime.
#' @importFrom  rgl translate3d extrude3d shade3d
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
                   st=NULL,alpha=1){
  # get length of STP_track
  n <- length(STP_track)
  # time_interval in seconds
  time_interval<- time_interval*60

  # if no zfacor provided calculate one based on longets spatial axis
  if(is.null(zfactor)){
    bbox<-bbox(STP_track)
    spatial_diff<-max(bbox[1,2]-bbox[1,1],bbox[2,2]-bbox[2,1])
    time_diff<-difftime(STP_track@endTime[[n]],STP_track@endTime[[1]],units = 'mins')
    zfac<-spatial_diff/as.numeric(time_diff)
  }else{
  zfac<-zfactor
}

  # if no st provided. start at first time of STP_track
  if (is.null(st)){
    st <- STP_track@endTime[[1]]
  }
  # time uncertainty in secs
  tu <- STP_track@rough_sets$time_uncertainty*60

  # loop through STPS
  for(i in 1:(length(STP_track)-1)){

  # create list of times for which PPAS need to be calculated
  times<-seq(STP_track@endTime[i],STP_track@endTime[i+1],(time_interval))

  # calculate PPAS
  suppressWarnings(PPAS<-lapply(times, function(x) {
    PPA(STP_track, x,quadsegs = 12)
  }))

  # create times and ppas tip and top of STP_track in case of time uncertainty
  if(tu>0 & (i==1 | i==(length(STP_track)-1))){
    STP1 <- STP_track[i:(i+1),'']
    STP1@rough_sets$time_uncertainty <- 0
    STP1@endTime[1] <- STP1@endTime[1] - tu
    STP1@endTime[2] <- STP1@endTime[2] + tu


    if(i==1){
      times1 <- seq(STP_track@endTime[1]-tu,STP_track@endTime[1]-time_interval,(time_interval))
    }else{
      times1 <- seq(STP_track@endTime[n]+time_interval,STP_track@endTime[n]+tu,(time_interval))
    }

    suppressWarnings(PPAS1<-lapply(times1, function(x) {
      calculate_PPA(STP1, x)
    }))
    if(i==1){
      PPAS<-c(PPAS1,PPAS)
      times<-c(times1,times)}
    else{
      PPAS<-c(PPAS,PPAS1)
      times<-c(times,times1)
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
  xx<-STP_coords$x
  yy<-STP_coords$y
  zz<-STP_coords$z

  # add orignal space-time points to STP. Only if no PPA could be calculated
  if(NAs[1]==F){
    xx<-c(STP_track@sp[i]@coords[1],xx)
    yy<-c(STP_track@sp[i]@coords[2],yy)
    zz<-c(tdif[1],zz)
  }# also if location uncertainty is 0 and time of space-time points is not in times
  if(tail(NAs,1)==F | STP_track@rough_sets$location_uncertainty==0){
    xx<-c(xx,STP_track@sp[(i+1)]@coords[1])
    yy<-c(yy,STP_track@sp[(i+1)]@coords[2])
    zz<-c(zz,(as.numeric(difftime(STP_track@endTime[(i+1)],st,units = 'min'))+tu/60)*zfac)
  }
  # matrix with all coordinates
  stp3d<-matrix(c(xx,yy,zz),ncol=3)

  at <- STP_track@connections$activity_time[i]*60
  # if activity time bigger than 0 take convexhull of stp parts separately
  if (at>0){
    # calculate halfway time
    tdiff<-abs(difftime(STP_track@endTime[i],STP_track@endTime[(i+1)],units = 'secs'))
    middle_time <- STP_track@endTime[i]+tdiff/2

    b1 <- as.numeric(difftime(middle_time-at/2,st,units='mins')*zfac)
    b2 <- as.numeric(difftime(middle_time+at/2,st,units='mins')*zfac)

    # split inot upper cone bottum cone and zone of no movement(part2)
    part1<-stp3d[stp3d[,3]<=(b1+time_interval/60*zfac),]
    part2<-stp3d[stp3d[,3]>=b1 & stp3d[,3]<=b2,]
    part3<-stp3d[stp3d[,3]>=(b2-time_interval/60*zfac),]

    # get convexhull
    conv1<-t(convhulln(part1))
    conv2<-t(convhulln(part2))
    conv3<-t(convhulln(part3))
    # plot the three parts
    rgl.triangles(part1[conv1,1],part1[conv1,2],part1[conv1,3],col=col,alpha=alpha)
    rgl.triangles(part2[conv2,1],part2[conv2,2],part2[conv2,3],col=col,alpha=alpha)
    rgl.triangles(part3[conv3,1],part3[conv3,2],part3[conv3,3],col=col,alpha=alpha)
  }else{
    # take convexhull and plot STP

    conv<-t(convhulln(stp3d))
    rgl.triangles(stp3d[conv,1],stp3d[conv,2],stp3d[conv,3],col=col,alpha=alpha)

  }
  }
  # return zfactor is it was not provided
  if(is.null(zfactor)){
return(zfac)}

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

