dir.create("~/Week10")
setwd("~/Week10/")
list.files(all.files = T)
objects()   # Should be empty.
rm(list=objects())
#

if (!require("pacman")) install.packages("pacman")
pacman::p_load(httr,EcoHydRology,curl,elevatr,raster,rgdal,
               data.table,foreign,maptools,dataRetrieval,gdistance, ggplot2)

make_usgs_gage_list=function(siteNo = "0205551460",
                             parameterCd = c("00060","00065"), # 60 discharge, 65 gage height https://waterdata.usgs.gov/nwis/uv?site_no=0205551460
                             # shorter time period to see individual waves
                             start.date = "2017-05-01",  # Not frozen to not frozen
                             end.date = "2017-11-01"){    # to still not frozen
  
  USGSlist=list()   # Organize the data in a nice list as in previous labs
  USGSlist[["flowdata"]]<- readNWISuv(siteNumbers = siteNo,parameterCd = parameterCd,startDate = start.date,endDate = end.date)
  head(USGSlist$flowdata)  # Note that we have 00060 and 00065...
  
  #
  # And of course we want to work in SI units so:
  USGSlist$flowdata$depth_m=USGSlist$flowdata$X_00065_00000*0.3048
  # m/ft depth
  USGSlist$flowdata$cms=USGSlist$flowdata$X_00060_00000*.02832
  # m3/ft3 flow
  #
  # Let's add in the USGS gage site information to the list and inspect
  USGSlist[["site"]]=readNWISsite(siteNo)
  head(USGSlist$site)
  class(USGSlist$site$dec_lat_va)
  #
  # Set the Manning Coefficient in the USGS Gage's Site Table
  #
  USGSlist$site$man_n=.035/1.49
  #
  # Create a SpatialPointsDataFrame out of the site dataframe in the USGS list
  coordinates(USGSlist$site)=~dec_long_va+dec_lat_va
  #
  
  return(USGSlist)
}

USGS02056000=make_usgs_gage_list(siteNo = "02056000")
USGS02055100=make_usgs_gage_list(siteNo = "02055100")
USGS02055000=make_usgs_gage_list(siteNo = "02055000")
USGS02054530=make_usgs_gage_list(siteNo = "02054530")

ab_ll=rbind(USGS02056000$site,
            USGS02055100$site,
            USGS02055000$site,
            USGS02054530$site)
class(ab_ll)
ab_ll@proj4string
proj4_utm = paste0("+proj=utm +zone=",
                   trunc((180+coordinates(USGS02055000$site)[1])/6+1), 
                   " +datum=WGS84 +units=m +no_defs")
print(proj4_utm)
# Lat/Lon (_ll) is much easier!
proj4_ll = "+proj=longlat"
crs_ll=CRS(proj4_ll)
crs_utm=CRS(proj4_utm)
proj4string(ab_ll)=proj4_ll
ab_utm=spTransform(ab_ll,crs_utm)
ab_utm@coords
mydem=get_aws_terrain(locations=ab_utm@coords, 
                      z = 12, prj = proj4_utm,expand=1)
#
# Lets plot the DEM and the gage locations so we can guess 
# what gages connect with what gages
#
plot(mydem)
plot(ab_utm,add=T)
text(ab_utm, labels=ab_utm@data$site_no, cex=0.6, font=2,pos=1)
#
# NOTE Since this is the upper most gage, we will calculate the distance from # this one to the next gage once we get the info for the next gage

# Set the starting and ending locations
# determine the river reach length and slope using the gdistance package.
#

# gdistance calculates the distance between given geometries

#1) 02055100 Tinker Creek at Daleville to 02056000 Roanoke River at Niagara
#2) 02054530 Roanoke River at Glenvar to 02055000 Roanoke River at Roanoke
#3) 02055000 Roanoke River at Roanoke to 02056000 Roanoke River at Niagara

Niagara=SpatialPoints(USGS02056000$site) # Down gradient site ROA River at Niagara
proj4string(Niagara)=proj4_ll
Niagara_utm=spTransform(Niagara,crs_utm)

##### Tinker Creek to Niagara #####
Tinker=SpatialPoints(USGS02055100$site)# Up gradient site Lick Run
proj4string(Tinker)=proj4_ll
Tinker_utm=spTransform(Tinker,crs_utm)
# Cut the DEM down to a more manageable size
cropmydem=crop(mydem,extend(extent(ab_utm),600))
cropmydem=trim(cropmydem)
cropmydem=cropmydem*1000.0
# Set up the weighting functions
altDiff <- function(x){x[2] - x[1]}
hd <- transition(cropmydem, altDiff, 8, symm=FALSE)
slope <- geoCorrection(hd)
adj <- adjacent(cropmydem, cells=1:ncell(cropmydem), pairs=TRUE, directions=8)
speed <- slope
speed[adj] <- 6 * exp(-3.5 * abs(slope[adj] + 0.05))
Conductance <- geoCorrection(speed)
# Find and plot the flow path
TinkertoNiagara <- shortestPath(Conductance, Tinker_utm, Niagara_utm, output="SpatialLines")
SpatialLinesLengths(TinkertoNiagara)

attach(USGS02055100)
site$L=SpatialLinesLengths(TinkertoNiagara) # km to m
site$L # reach length in m
#
#
# Getting slope, we will extract the slope for points A and B from the DEM and # divide the difference by the length in m, this gives us a much better 
# estimate of slope than taking the point slopes at the gage site
#
site$slope=(extract(mydem,Tinker_utm)-
                             extract(mydem,Niagara_utm))/site$L
site$slope

flowdata$Niagara=(site$man_n*
                             flowdata$cms)/(flowdata$depth_m^(5/3)*
                                                             sqrt(site$slope))
head(flowdata)

# Lets look at how B changes with flow.    
#plot(flowdata$dateTime,flowdata$Niagara, main="LICK RUN TO ROANOKE RIVER AT NIAGARA, VA")
# Does this seem reasonable (...like order of magnitude reasonable)? You can 
# perform a quick and dirty check using google earth and measuring the channel 
# width in a few places.
#
#plot(flowdata$cms,flowdata$depth_m, main="LICK RUN TO ROANOKE RIVER AT NIAGARA, VA")

flowdata$ck =
  5/3*sqrt(site$slope)/site$man_n*
  (flowdata$depth_m^(2/3))
#
flowdata$dt =
  site$L/flowdata$ck

plot(flowdata$dateTime,flowdata$dt)
# time for the wave to move down the channel

flowdata$outTime=flowdata$dateTime+
  flowdata$dt

# Find beginning of  Waves
flowdata$newwave=
  flowdata$cms *1.1 <
  data.table::shift(flowdata$cms)
# 10% larger
summary(flowdata$newwave)
# when it switches to true, there is an increase in flow (wave)
# Add plot of the point found
len=length(flowdata$newwave)
flowdata$newwave[is.na(flowdata$newwave)]=F
# Removes repeated finds by going through loop backwords
for (i in seq(len,2)){
  print(i)
  if(flowdata$newwave[i]==T &
     flowdata$newwave[i-1]==T){
    flowdata$newwave[i]=F
  }
}
plot(flowdata$dateTime,flowdata$cms,type="l")
points(flowdata$dateTime[flowdata$newwave],
       flowdata$cms[flowdata$newwave],col=2)

# Find the time locations where waves begin
which(flowdata$newwave == TRUE)

ggplot() +
  geom_line(data=flowdata, aes(x=dateTime, y=cms, colour="Start Time")) +
  geom_line(data=flowdata, aes(x=outTime, y=cms, colour = "Out Time")) +
  xlab("Date") + ylab("Flow (cms)") + 
  xlim(c(flowdata$dateTime[1109],flowdata$dateTime[1109+200])) +
  ylim(c(0,3)) + ggtitle("Tinker Creek at Daleville to ROA River at Niagara")

detach(USGS02055100)

##### Glenvar to Niagara #####

Glenvar=SpatialPoints(USGS02054530$site)# Up gradient site Lick Run
proj4string(Glenvar)=proj4_ll
Glenvar_utm=spTransform(Glenvar,crs_utm)
# Cut the DEM down to a more manageable size
cropmydem=crop(mydem,extend(extent(ab_utm),600))
cropmydem=trim(cropmydem)
cropmydem=cropmydem*1000.0
# Set up the weighting functions
altDiff <- function(x){x[2] - x[1]}
hd <- transition(cropmydem, altDiff, 8, symm=FALSE)
slope <- geoCorrection(hd)
adj <- adjacent(cropmydem, cells=1:ncell(cropmydem), pairs=TRUE, directions=8)
speed <- slope
speed[adj] <- 6 * exp(-3.5 * abs(slope[adj] + 0.05))
Conductance <- geoCorrection(speed)
# Find and plot the flow path
GlenvartoNiagara <- shortestPath(Conductance, Glenvar_utm, Niagara_utm, output="SpatialLines")
SpatialLinesLengths(GlenvartoNiagara)

attach(USGS02054530)
site$L=SpatialLinesLengths(GlenvartoNiagara) # km to m
site$L # reach length in m
#
#
# Getting slope, we will extract the slope for points A and B from the DEM and # divide the difference by the length in m, this gives us a much better 
# estimate of slope than taking the point slopes at the gage site
#
site$slope=(extract(mydem,Glenvar_utm)-
              extract(mydem,Niagara_utm))/site$L
site$slope

flowdata$Niagara=(site$man_n*
                    flowdata$cms)/(flowdata$depth_m^(5/3)*
                                     sqrt(site$slope))
head(flowdata)

# Lets look at how B changes with flow.    
#plot(flowdata$dateTime,flowdata$Niagara, main="LICK RUN TO ROANOKE RIVER AT NIAGARA, VA")
# Does this seem reasonable (...like order of magnitude reasonable)? You can 
# perform a quick and dirty check using google earth and measuring the channel 
# width in a few places.
#
#plot(flowdata$cms,flowdata$depth_m, main="LICK RUN TO ROANOKE RIVER AT NIAGARA, VA")

flowdata$ck =
  5/3*sqrt(site$slope)/site$man_n*
  (flowdata$depth_m^(2/3))
#
flowdata$dt =
  site$L/flowdata$ck

plot(flowdata$dateTime,flowdata$dt)
# time for the wave to move down the channel

flowdata$outTime=flowdata$dateTime+
  flowdata$dt

# Find beginning of  Waves
flowdata$newwave=
  flowdata$cms *1.1 <
  data.table::shift(flowdata$cms)
# 10% larger
summary(flowdata$newwave)
# when it switches to true, there is an increase in flow (wave)
# Add plot of the point found
len=length(flowdata$newwave)
flowdata$newwave[is.na(flowdata$newwave)]=F
# Removes repeated finds by going through loop backwords
for (i in seq(len,2)){
  print(i)
  if(flowdata$newwave[i]==T &
     flowdata$newwave[i-1]==T){
    flowdata$newwave[i]=F
  }
}
plot(flowdata$dateTime,flowdata$cms,type="l")
points(flowdata$dateTime[flowdata$newwave],
       flowdata$cms[flowdata$newwave],col=2)

# Find the time locations where waves begin
which(flowdata$newwave == TRUE)

ggplot() +
  geom_line(data=flowdata, aes(x=dateTime, y=cms, colour="Start Time")) +
  geom_line(data=flowdata, aes(x=outTime, y=cms, colour = "Out Time")) +
  xlab("Date") + ylab("Flow (cms)") + 
  xlim(c(flowdata$dateTime[1109],flowdata$dateTime[1109+200])) +
  ylim(c(0,3)) + ggtitle("ROA River at Glenvar to ROA River at Niagara")

detach(USGS02054530)

##### Roanoke to Niagara #####

Roanoke=SpatialPoints(USGS02055000$site)# Up gradient site Lick Run
proj4string(Roanoke)=proj4_ll
Roanoke_utm=spTransform(Roanoke,crs_utm)
# Cut the DEM down to a more manageable size
cropmydem=crop(mydem,extend(extent(ab_utm),600))
cropmydem=trim(cropmydem)
cropmydem=cropmydem*1000.0
# Set up the weighting functions
altDiff <- function(x){x[2] - x[1]}
hd <- transition(cropmydem, altDiff, 8, symm=FALSE)
slope <- geoCorrection(hd)
adj <- adjacent(cropmydem, cells=1:ncell(cropmydem), pairs=TRUE, directions=8)
speed <- slope
speed[adj] <- 6 * exp(-3.5 * abs(slope[adj] + 0.05))
Conductance <- geoCorrection(speed)
# Find and plot the flow path
RoanoketoNiagara <- shortestPath(Conductance, Roanoke_utm, Niagara_utm, output="SpatialLines")
SpatialLinesLengths(RoanoketoNiagara)

attach(USGS02055000)
site$L=SpatialLinesLengths(RoanoketoNiagara) # km to m
site$L # reach length in m
#
#
# Getting slope, we will extract the slope for points A and B from the DEM and # divide the difference by the length in m, this gives us a much better 
# estimate of slope than taking the point slopes at the gage site
#
site$slope=(extract(mydem,Roanoke_utm)-
              extract(mydem,Niagara_utm))/site$L
site$slope

flowdata$Niagara=(site$man_n*
                    flowdata$cms)/(flowdata$depth_m^(5/3)*
                                     sqrt(site$slope))
head(flowdata)

# Lets look at how B changes with flow.    
#plot(flowdata$dateTime,flowdata$Niagara, main="LICK RUN TO ROANOKE RIVER AT NIAGARA, VA")
# Does this seem reasonable (...like order of magnitude reasonable)? You can 
# perform a quick and dirty check using google earth and measuring the channel 
# width in a few places.
#
#plot(flowdata$cms,flowdata$depth_m, main="LICK RUN TO ROANOKE RIVER AT NIAGARA, VA")

flowdata$ck =
  5/3*sqrt(site$slope)/site$man_n*
  (flowdata$depth_m^(2/3))
#
flowdata$dt =
  site$L/flowdata$ck

plot(flowdata$dateTime,flowdata$dt)
# time for the wave to move down the channel

flowdata$outTime=flowdata$dateTime+
  flowdata$dt

# Find beginning of  Waves
flowdata$newwave=
  flowdata$cms *1.1 <
  data.table::shift(flowdata$cms)
# 10% larger
summary(flowdata$newwave)
# when it switches to true, there is an increase in flow (wave)
# Add plot of the point found
len=length(flowdata$newwave)
flowdata$newwave[is.na(flowdata$newwave)]=F
# Removes repeated finds by going through loop backwords
for (i in seq(len,2)){
  print(i)
  if(flowdata$newwave[i]==T &
     flowdata$newwave[i-1]==T){
    flowdata$newwave[i]=F
  }
}
plot(flowdata$dateTime,flowdata$cms,type="l")
points(flowdata$dateTime[flowdata$newwave],
       flowdata$cms[flowdata$newwave],col=2)

# Find the time locations where waves begin
which(flowdata$newwave == TRUE)

ggplot() +
  geom_line(data=flowdata, aes(x=dateTime, y=cms, colour="Start Time")) +
  geom_line(data=flowdata, aes(x=outTime, y=cms, colour = "Out Time")) +
  xlab("Date") + ylab("Flow (cms)") + 
  xlim(c(flowdata$dateTime[1109],flowdata$dateTime[1109+200])) +
  ylim(c(0,3)) + ggtitle("ROA River at Roanoke to ROA River at Niagara")

detach(USGS02055000)
