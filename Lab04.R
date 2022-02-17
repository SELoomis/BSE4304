### BSE 4304 - Lab 04
### Sarah Loomis
### 18 February 2022

## Cleaning up and setting up directory
options(repos ="http://cran.us.r-project.org")  # required to get latest libs
objects()
if (!require("pacman")) install.packages("pacman")
pacman::p_load(elevatr,soilDB,rgdal,raster)
pacman::p_load(EcoHydRology,rnoaa,curl,httr,ggplot2)
rm(list=objects())
setwd("~")
dir.create("~/Lab04")
setwd("~/Lab04/")

## Source code for soil wetting/drying/etc functions
url="https://raw.githubusercontent.com/vtdrfuka/BSE5304_2022/main/Lab04.R"
download.file(url,"Lab04.R")
browseURL("https://raw.githubusercontent.com/vtdrfuka/BSE5304_2022/main/functions/TMWBmodel.R")
source("https://raw.githubusercontent.com/vtdrfuka/BSE5304_2022/main/functions/TMWBmodel.R")

## Downloading soils data from Web Soil Survey
url="https://websoilsurvey.sc.egov.usda.gov/DSD/Download/AOI/tjhr1b3f3qjq4lcfn3ludyx0/wss_aoi_2022-02-17_09-39-44.zip"
download.file(url,"mysoil.zip")
unzip("mysoil.zip")

## Getting data from USGS Gage
myflowgage_id="01415460"
myflowgage=get_usgs_gage(myflowgage_id,begin_date = "2015-01-01",end_date = "2022-03-01")
myflowgage$flowdata$Qmm = myflowgage$flowdata$flow/myflowgage$area/10^3
  # change data from m3/day to mm/day

## Work with soils data
mysoil=readOGR("wss_aoi_2022-02-17_09-39-44/spatial/soilmu_a_aoi.shp")    
  # Explore the mysoil dataset which is returned
mybbox=c(mysoil@bbox)
  # First associate mukey with cokey from component
mysoil$mukey=mysoil$MUKEY  # or rename the column
mukey_statement = format_SQL_in_statement(unique(mysoil$mukey))
q_mu2co = paste("SELECT mukey,cokey FROM component WHERE mukey IN ", mukey_statement, sep="")
mu2co = SDA_query(q_mu2co)
  # Second associate cokey with ksat_r,awc_r,hzdepb_r from chorizon
cokey_statement = format_SQL_in_statement(unique(mu2co$cokey))
q_co2ch = paste("SELECT cokey,ksat_r,awc_r,hzdepb_r  FROM chorizon WHERE cokey IN ", cokey_statement, sep="")
co2ch = SDA_query(q_co2ch)
  # Last, bring them back together, and aggregate based on max values of ksat_r,awc_r, and hzdepb_r
mu2ch=merge(mu2co,co2ch)
mu2ch$awc_r[is.na(mu2ch$awc_r)]=0
summary(mu2ch)
  # should be no NAs
mu2chmax=aggregate(mu2ch,list(mu2ch$mukey),max)
mu2chmax$depthmm=mu2chmax$hzdepb_r*10
  # change depth to restrictive layer from cm to mm
mu2chmax$AWC=mu2chmax$awc_r*mu2chmax$depthmm
  # calculate AWC (depth in mm) from percentage to depth by multiplying by depth to restrictive layer (mm)

## Elevation Data
proj4_ll = "+proj=longlat"
proj4string(mysoil) = proj4_ll
mydem=get_elev_raster(locations=mysoil, 
                      z = 11, prj =proj4string(mysoil) ,
                      src ="aws",clip="bbox",expand = 0.001)

summary(terrain(mydem, opt='slope',unit = "degrees"))
plot(terrain(mydem,opt='slope',unit="degrees"))
lines(mysoil,col="black")
  # plot the slope with the soil layer for the drainage area

## Weather Data
stns=meteo_distance(
  station_data=ghcnd_stations(),
  lat=myflowgage$declat,
  long=myflowgage$declon,
  units = "deg",
  radius = 30,
  limit = NULL
)
 
WXStn=stns[stns$element=="TMAX"&stns$last_year>=2021,]$id[1]
WXData=meteo_pull_monitors(
  monitors=WXStn,
  keep_flags = FALSE,
  date_min = "2016-01-01",
  date_max = NULL,
  var = c("TMAX","TMIN","PRCP") 
)
summary(WXData)  

## Create an aligned modeldata data frame to build our model in
modeldata=merge(WXData,myflowgage$flowdata,by.x="date",by.y="mdate")
modeldata$MaxTemp=modeldata$tmax/10 # Converting to C
modeldata$MinTemp=modeldata$tmin/10 # Converting to C
modeldata$P=modeldata$prcp/10 # Converting to mm
modeldata$P[is.na(modeldata$P)]=0
modeldata$MinTemp[is.na(modeldata$MinTemp)]=0
modeldata$MaxTemp[is.na(modeldata$MaxTemp)]= modeldata$MinTemp[is.na(modeldata$MaxTemp)] +1
modeldata$MaxTemp[modeldata$MaxTemp<=modeldata$MinTemp]= modeldata$MinTemp[modeldata$MaxTemp<=modeldata$MinTemp]+1
modeldata$AvgTemp=(modeldata$MaxTemp+modeldata$MaxTemp)/2.0
summary(modeldata)
  # should be no NAs
TMWB=modeldata

## Calibrating the parameters one at a time
for (fcres in seq(.1,.5,.1)){
  TMWBnew=TMWBmodel(TMWB=TMWB,fcres=fcres)
  print(paste(fcres,NSE(TMWBnew$Qmm,TMWBnew$Qpred)))
}
  # fcres=0.1 gives highest NSE for this location
for (SFTmp in seq(-5,20)){
  TMWBnew=TMWBmodel(TMWB=TMWB,fcres=.1,SFTmp = SFTmp)
  print(paste(SFTmp,NSE(TMWBnew$Qmm,TMWBnew$Qpred)))
}
  # SFTmp=-5 gives highest NSE for this location
for (bmlt6 in seq(1.4,7,0.2)){
  TMWBnew=TMWBmodel(TMWB=TMWB,fcres=.1,SFTmp = -5,bmlt6=bmlt6)
  print(paste(bmlt6,NSE(TMWBnew$Qmm,TMWBnew$Qpred)))
}
  # bmlt6 gave the same NSE values - picked 2.5
for (bmlt12 in seq(1.4,7,0.2)){
  TMWBnew=TMWBmodel(TMWB=TMWB,fcres=.1,SFTmp = -5,bmlt6=2.5,bmlt12=bmlt12)
  print(paste(bmlt12,NSE(TMWBnew$Qmm,TMWBnew$Qpred)))
}
  # bmlt12 gave the same NSE values - picked 2.5
for (Tlag in seq(0,1,0.1)){
  TMWBnew=TMWBmodel(TMWB=TMWB,fcres=.1,SFTmp = -5,bmlt6=2.5,bmlt12=2.5,Tlag=Tlag)
  print(paste(Tlag,NSE(TMWBnew$Qmm,TMWBnew$Qpred)))
}
  # Tlag provided the same NSE for all values - choose 0.5
# Based on calculated AWC data, estimate AWC values to be between 135.2 and 732
for(AWCval in seq(135,732,50)){
  TMWBnew=TMWBmodel(TMWB=TMWB,fcres=.1,SFTmp = -5,bmlt6=2.5,bmlt12=2.5,Tlag = .5,AWCval=AWCval)
  print(paste(AWCval,NSE(TMWBnew$Qmm,TMWBnew$Qpred)))
}
  # AWCval=135 gives highest NSE for this location
  # Best result for "Terry Clove Kill near De Lancey, NY" NSE = 0.185
TMWBnew=TMWBmodel(TMWB=TMWB,fcres=.1,SFTmp = -5,bmlt6=2.5,bmlt12=2.5,Tlag = .5,AWCval = 135.2,Slope=0)

# For a simple spatial model, use TMWB to initialize 3 slope components for Top, Mid, and Bottom of the hillside. 
# Based on looking at the plot of slope and soils data, estimate a hillslope to have topslope=15,midslope=25,bottomslope=5 
TopSlope=TMWBmodel(TMWB=TMWBnew,AWCval = 200,Slope = 15)

MidSlope=TopSlope
MidSlope$P=MidSlope$P+MidSlope$Excess
MidSlope=TMWBmodel(TMWB=MidSlope,AWCval = 600,Slope = 25)

BotSlope=MidSlope
BotSlope$P=BotSlope$P+BotSlope$Excess
BotSlope=TMWBmodel(TMWB=BotSlope,AWCval = 400,Slope = 5)

## Plot for AW
ggplot() +
  geom_line(data=BotSlope,aes(x=date, y=AW, colour ="BotSlope, Slope=5%"), size=1) +
  geom_line(data=MidSlope,aes(x=date, y=AW, colour ="MidSlope, Slope=25%"), size=1) +
  geom_line(data=TopSlope,aes(x=date, y=AW, colour ="TopSlope, Slope=15%"), size=1) +
  xlab("Date") + ylab("AW (mm)") + ggtitle("USGS Gage 01415460 – Terry Clove Kill near De Lancey, NY")

## Plot for Excess
ggplot() +
  geom_line(data=BotSlope,aes(x=date, y=Excess, colour ="BotSlope, Slope=5%"), size=1) +
  geom_line(data=MidSlope,aes(x=date, y=Excess, colour ="MidSlope, Slope=25%"), size=1) +
  geom_line(data=TopSlope,aes(x=date, y=Excess, colour ="TopSlope, Slope=15%"), size=1) +
  xlab("Date") + ylab("Excess (mm)") + ggtitle("USGS Gage 01415460 – Terry Clove Kill near De Lancey, NY")
