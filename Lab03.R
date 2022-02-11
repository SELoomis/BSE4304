### BSE 4304 Lab 03
### Author: Sarah Loomis
### Date: February 11th 2022

## Loading in Libraries
if (!require("pacman")) install.packages("pacman")
pacman::p_load(rnoaa,EcoHydRology,lattice,ggplot2)

## Calculating NSE, set up function for calibration
NSE=function(Yobs,Ysim){
  return(1-sum((Yobs-Ysim)^2, na.rm=TRUE)/sum((Yobs-mean(Yobs, na.rm=TRUE))^2, na.rm=TRUE))
}

## USGS Gage and Weather Data
myflowgage_id="04235000"
myflowgage=get_usgs_gage(myflowgage_id,
                         begin_date="2016-01-01",end_date="2022-02-09")
myflowgage$flowdata$flowmm=myflowgage$flowdata$flow/myflowgage$area/10^3
# normalize flow data to be in mm of flow

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
WXData$prcp=WXData$prcp/10
# adjusts tenths of mm to mm precipitation

## Merging and creating TMWB dataset (same as Lab02)
modeldata=merge(WXData, myflowgage$flowdata, by.x="date", by.y="mdate")
names(modeldata)[names(modeldata) == "flowmm"] <- "Qmm"
names(modeldata)[names(modeldata) == "prcp"] <- "P"
names(modeldata)[names(modeldata) == "tmax"] <- "MaxTemp"
names(modeldata)[names(modeldata) == "tmin"] <- "MinTemp"
modeldata$MaxTemp=modeldata$MaxTemp/10
modeldata$MinTemp=modeldata$MinTemp/10
modeldata$MinTemp[is.na(modeldata$MinTemp)]=0
modeldata$MaxTemp[is.na(modeldata$MaxTemp)]=modeldata$MinTemp[is.na(modeldata$MaxTemp)] +1
modeldata$P[is.na(modeldata$P)]=0
# adjusts tenths of degC to degC
TMWB=modeldata

## Using SnowMelt function to test different slopes/aspects
attach(TMWB)
SNO_Energy1=SnowMelt(date, P, MaxTemp-3, MinTemp-3, myflowgage$declat, 
                    slope = atan(0/100),
                    aspect = 0*(pi/180), tempHt = 1, windHt = 2, groundAlbedo = 0.25,
                    SurfEmissiv = 0.95, windSp = 2, forest = 0, startingSnowDepth_m = 0,
                    startingSnowDensity_kg_m3=450)
SNO_Energy2=SnowMelt(date, P, MaxTemp-3, MinTemp-3, myflowgage$declat, 
                     slope = atan(10/100),
                     aspect = 0*(pi/180), tempHt = 1, windHt = 2, groundAlbedo = 0.25,
                     SurfEmissiv = 0.95, windSp = 2, forest = 0, startingSnowDepth_m = 0,
                     startingSnowDensity_kg_m3=450)
SNO_Energy3=SnowMelt(date, P, MaxTemp-3, MinTemp-3, myflowgage$declat, 
                     slope = atan(10/100),
                     aspect = 180*(pi/180), tempHt = 1, windHt = 2, groundAlbedo = 0.25,
                     SurfEmissiv = 0.95, windSp = 2, forest = 0, startingSnowDepth_m = 0,
                     startingSnowDensity_kg_m3=450)
SNO_Energy4=SnowMelt(date, P, MaxTemp-3, MinTemp-3, myflowgage$declat, 
                     slope = atan(45/100),
                     aspect = 315*(pi/180), tempHt = 1, windHt = 2, groundAlbedo = 0.25,
                     SurfEmissiv = 0.95, windSp = 2, forest = 0, startingSnowDepth_m = 0,
                     startingSnowDensity_kg_m3=450)
SNO_Energy5=SnowMelt(date, P, MaxTemp-3, MinTemp-3, myflowgage$declat, 
                     slope = atan(45/100),
                     aspect = 225*(pi/180), tempHt = 1, windHt = 2, groundAlbedo = 0.25,
                     SurfEmissiv = 0.95, windSp = 2, forest = 0, startingSnowDepth_m = 0,
                     startingSnowDensity_kg_m3=450)
detach(TMWB)

## Plotting for HW 1
ggplot() +
  geom_line(data=SNO_Energy1,aes(x=Date, y=SnowWaterEq_mm, colour ="0% Slope, N Facing"), size=1) +
  geom_line(data=SNO_Energy2,aes(x=Date, y=SnowWaterEq_mm, colour ="10% Slope, N Facing"), size=1) +
  geom_line(data=SNO_Energy3,aes(x=Date, y=SnowWaterEq_mm, colour ="10% Slope, S Facing"), size=1) +
  geom_line(data=SNO_Energy4,aes(x=Date, y=SnowWaterEq_mm, colour ="45% Slope, NW Facing"), size=1) +
  geom_line(data=SNO_Energy5,aes(x=Date, y=SnowWaterEq_mm, colour ="45% Slope, SW Facing"), size=1) +
  ylab("Snow Water Equivalent (mm)") + ggtitle("USGS Gage 04235000 – Canandaigua Outlet at Chapin NY")
  
## Define function for HW 2 calibration
SAF=function(fcres=0.3,SFTmp = 1,bmlt6 = 4.5,bmlt12 = 0.0,Tlag = 1){
  
## Calculate melt factor
#SFTmp = 1 #snowfall temp in degC
#bmlt6 = 4.5 #melt factor for snow on 6/21 (mm water / degC*day)
#bmlt12 = 0  #melt factor for snow on 12/21 (mm water / degC*day)
#Tlag = 1  #snow pack temperature lag factor
Tmlt = SFTmp  # Assumed to be same as SnowFall Temperature
TMWB$AvgTemp=(TMWB$MaxTemp + TMWB$MinTemp)/2
TMWB$bmlt = (bmlt6 + bmlt12)/2 + (bmlt6 - bmlt12)/2 *  sin(2*pi/365*(julian(TMWB$date,origin = as.Date("2000-01-01"))-81))
  # calculates melt factor based on eqn. 1:2.5.3

## Initialize SNO, Tsno as well as the first values of each
TMWB$SNO = 0  # Snow Depth (mm)
TMWB$Tsno = 0  # Snow Temp (C)
TMWB$SNOmlt = 0  # Snow Melt (mm)
attach(TMWB)
for (t in 2:length(date)){
  Tsno[t]= Tsno[t-1] * (1.0-Tlag) +  AvgTemp[t] * Tlag
    # calculates snow pack temp based on eqn 1:2.5.1
  if(AvgTemp[t] < SFTmp){
    SNO[t]= SNO[t-1] + P[t]
      # if the average temp is less than snowfall temp, snow depth is previous day's snow depth plus current precipitation
  }  else {
      # if the average temp is greater than snowfall temp
    SNOmlt[t]= bmlt[t] * SNO[t-1] * ((Tsno[t]+MaxTemp[t])/2 - Tmlt) 
      # calculates snowmelt based on eqn 1:2.5.2
    SNOmlt[t]= min(SNOmlt[t],SNO[t-1])
      # snowmelt is the smaller of snowmelt or previous day's snow depth
    SNO[t]= SNO[t-1] -SNOmlt[t]
      # snow depth = previous day's snow depth minus today's snowmelt
  }
}
detach(TMWB)
TMWB$Tsno=Tsno
TMWB$SNO=SNO
TMWB$SNOmlt=SNOmlt
rm(list=c("SNO", "SNOmlt", "Tsno"))

## Estimate PET using Albedo
TMWB$Albedo=.23
  # albedo = 0.23 if there's no snow
TMWB$Albedo[TMWB$SNO>0]=.95
  # albedo = 0.95 if there is snow
attach(TMWB)
PET=PET_fromTemp(Jday=(1+as.POSIXlt(date)$yday),Tmax_C = MaxTemp,Tmin_C = MinTemp,albedo=Albedo,lat_radians = myflowgage$declat*pi/180) * 1000
TMWB$PET=PET
plot(date,PET)
detach(TMWB)
rm(list=c("PET"))

## Initializaing AWC, Net Precip, ET, AW, and Excess
#myflowgage$FldCap=.45
#myflowgage$WiltPt=.15
#myflowgage$Z=1000
TMWB$AWC=350 
TMWB$dP = 0 
TMWB$ET = 0 
TMWB$AW = 0 
TMWB$Excess = 0

soilwetting<-function(AWprev,dP_func,AWC_func){
  AW_func<-AWprev+dP_func
  excess_func<-0.0
  c(AW_func,excess_func)
} 

soildrying<-function(AWprev,dP_func,AWC_func){
  AW_func=AWprev*exp(dP_func/AWC_func)
  excess_func<-0.0
  c(AW_func,excess_func)
}
# soil_wetting_above_capacity function
soil_wetting_above_capacity<-function(AWprev,dP_func,AWC_func){
  AW_func<-AWC_func
  excess_func<-AWprev+dP_func-AWC_func
  c(AW_func,excess_func)
}

## Loop to calculate AW and Excess, this time using PET calculated earlier
attach(TMWB)
for (t in 2:length(AW)){
  ET[t] = min (AW[t-1],PET[t])
  ET[t] = (AW[t-1]/AWC[t-1])*PET[t] 
  if(AvgTemp[t] >= SFTmp){
    dP[t] = P[t] - ET[t] + SNOmlt[t] 
  }  else {
    dP[t] = ET[t]
  }
  if (dP[t]<=0) {
    values<-soildrying(AW[t-1],dP[t],AWC[t])
  } else if((dP[t]>0) & (AW[t-1]+dP[t])<=AWC[t]) {
    values<-soilwetting(AW[t-1],dP[t],AWC[t])
  } else {
    values<-soil_wetting_above_capacity(AW[t-1],dP[t],AWC[t])
  }
  AW[t]<-values[1]
  Excess[t]<-values[2]
}
TMWB$AW=AW
TMWB$Excess=Excess
TMWB$dP=dP
rm(list=c("AW","dP","ET", "Excess"))
detach(TMWB)

## Calculating predicted flow and storage
TMWB$Qpred=NA
TMWB$Qpred[1]=0
TMWB$S=NA
TMWB$S[1]=0
attach(TMWB)
#fcres=.3
for (t in 2:length(date)){
  S[t]=S[t-1]+Excess[t]     
  Qpred[t]=fcres*S[t]
  S[t]=S[t]-Qpred[t]
}
TMWB$S=S
TMWB$Qpred=Qpred 
detach(TMWB) 
rm(list=c("Qpred","S"))

return(NSE(TMWB$Qmm,TMWB$Qpred))
}
