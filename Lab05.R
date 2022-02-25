### BSE 4304 - Lab 05
### Sarah Loomis
### 25 February 2022 

## Cleaning up and setting up directory
options(repos ="http://cran.us.r-project.org")  # required to get latest libs
objects()
if (!require("pacman")) install.packages("pacman")
pacman::p_load(elevatr,soilDB,rgdal,raster)
pacman::p_load(EcoHydRology,rnoaa,curl,httr,ggplot2)
rm(list=objects())
#dir.create("~/Lab05")
setwd("~/Lab05/")
source("https://raw.githubusercontent.com/vtdrfuka/BSE5304_2022/main/functions/TMWBmodel.R")
source("https://raw.githubusercontent.com/SELoomis/BSE4304/main/Lab05_TMWB.R")
  # this sources the code for the previous Lab from GitHub to have the TMWB data to compare with CN data

## NSE Function
NSE=function(Yobs,Ysim){
  return(1-sum((Yobs-Ysim)^2, na.rm=TRUE)/sum((Yobs-mean(Yobs, na.rm=TRUE))^2, na.rm=TRUE))
}

## Curve Number Function
CNmodel<-function(CNmodeldf, CNavg = 75,IaFrac = 0.05,fnc_slope=0, 
                  fnc_aspect=0,func_DAWC=.3,func_z=1000,fnc_fcres=.3) {
  
  # Energy Balance based Snow Accumulation 
  # and Melt model from the EcoHydRology package.
  attach(CNmodeldf)
  SNO_Energy=SnowMelt(date, P, MaxTemp-3, MinTemp-3, myflowgage$declat, 
                      slope = fnc_slope, aspect = fnc_aspect, tempHt = 1, 
                      windHt = 2, groundAlbedo = 0.25,SurfEmissiv = 0.95, windSp = 2, 
                      forest = 0, startingSnowDepth_m = 0,startingSnowDensity_kg_m3=450)
  # We will update the -3 in the above to be a lapse rate adjustment
  detach(CNmodeldf)
  CNmodeldf$SNO=SNO_Energy$SnowWaterEq_mm
  CNmodeldf$SNOmlt=SNO_Energy$SnowMelt_mm
  CNmodeldf$SnowfallWatEq_mm=SNO_Energy$SnowfallWatEq_mm
  CNmodeldf$SnowMelt_mm=SNO_Energy$SnowMelt_mm
  attach(CNmodeldf)
  CNmodeldf$Albedo=.23
  CNmodeldf$Albedo[CNmodeldf$SNO>0]=.95
  PET=PET_fromTemp(Jday=(1+as.POSIXlt(date)$yday),
                   Tmax_C = MaxTemp,Tmin_C = MinTemp,
                   lat_radians = myflowgage$declat*pi/180) * 1000
  CNmodeldf$PET=PET
  detach(CNmodeldf)
  rm(list="PET")
  
  CNmodeldf$AWC=func_DAWC*func_z
  # Oh, this we want to vary some of these around our watershed!
  CNmodeldf$dP = 0 # Initializing Net Precipitation
  CNmodeldf$ET = 0 # Initializing ET
  CNmodeldf$AW = 0 # Initializing AW
  CNmodeldf$Excess = 0 # Initializing Excess
  CNmodeldf$S =0 # Initializing S
  CNmodeldf$Qpred=0 # Initializing Qpred
  attach(CNmodeldf)
  SSCNavg=(1000/CNavg-10)*25.4
  SSCN=SoilStorage(S_avg=SSCNavg, field_capacity=func_DAWC*.9,
                   soil_water_content=0.1*func_DAWC, porosity=func_DAWC)
  Ia_init=IaFrac*SSCN   
  CNmodeldf$CNavg = CNavg
  CNmodeldf$SSCNavg = SSCNavg
  CNmodeldf$SSCN = SSCN
  detach(CNmodeldf)
  rm(list=c("CNavg", "SSCN", "SSCNavg"))
  CNmodeldf$Ia = Ia_init
  attach(CNmodeldf)
  # Those processes that are dependant on prior days conditions, we run as a 
  # loop through each of the days.
  for (t in 2:length(AW)){
    ET[t] = AW[t-1]/AWC[t-1]*PET[t]
    # Calculating Net Precipitation which adds in slope above's Excess
    dP[t] = SNO_Energy$Rain_mm[t] - ET[t] + 
      SNO_Energy$SnowMelt_mm[t]    # CN Solution
    # Is the soil saturated, and thus can't take more dP? 
    if (AW[t-1] + dP[t]>=AWC[t]){
      Excess[t]=AW[t-1] + dP[t] -AWC[t]
      AW[t]=AWC[t]
      # Otherwise, if dP is less than the initial abstraction? 
      # https://en.wikipedia.org/wiki/Runoff_curve_number#Definition
    } else if (dP[t]<=Ia[t]) {
      Excess[t]=0.0
      AW[t]=AW[t-1] + dP[t]
    } else {
      Excess[t]=(dP[t]-Ia[t])^2/(dP[t]-Ia[t]+SSCN[t])
      AW[t]=AW[t-1] + dP[t] -Excess[t]
    }
    S[t]=S[t-1]+Excess[t]
    Qpred[t]=fnc_fcres*S[t]
    S[t]=S[t]-Qpred[t]
  }
  CNmodeldf$ET=ET
  CNmodeldf$dP=dP
  CNmodeldf$AW=AW
  CNmodeldf$Excess=Excess
  CNmodeldf$S=S
  CNmodeldf$Qpred=Qpred # UPDATE vector BEFORE DETACHING
  rm(list=c("AW", "dP", "ET", "Excess", "Qpred", "S"))
  detach(CNmodeldf)
  return(CNmodeldf)
}

## Curve Number calculations for hillslopes
TopSlopeCN=modeldata
MidSlopeCN=modeldata
BotSlopeCN=modeldata
  # Call the new CNmodel() function with Top,Mid,BotSlope HRU objects, 
  # passing the Qpred into the lower HRUs HillslopeAboveExcess (as area scaled flow)
TopSlopeCN = CNmodel(CNmodeldf = TopSlopeCN, CNavg = 60,IaFrac=0.05,fnc_slope=0,
                     fnc_aspect=0,func_DAWC=.3,
                     func_z=500,fnc_fcres=.3)
MidSlopeCN$P=TopSlopeCN$Excess+MidSlopeCN$P
# Higher slope, medium ksat, fcres=0.5 
MidSlopeCN = CNmodel(CNmodeldf = MidSlopeCN, CNavg = 60,IaFrac=0.05,fnc_slope=0, 
                     fnc_aspect=0,func_DAWC=.3,
                     func_z=750,fnc_fcres=.5)
# Low Slope and lowest ksat, fcres=0.2
BotSlopeCN$P=MidSlopeCN$Excess+BotSlopeCN$P
BotSlopeCN = CNmodel(CNmodeldf = BotSlopeCN, CNavg = 60,IaFrac=0.05,fnc_slope=0, 
                     fnc_aspect=0,func_DAWC=.3,
                     func_z=1000,fnc_fcres=.2)

## Plot TMWB Storage
p1 <- ggplot() +
  geom_line(data=BotSlope,aes(x=date, y=S, colour ="BotSlope"), size=1) +
  geom_line(data=MidSlope,aes(x=date, y=S, colour ="TopSlope"), size=1) +
  geom_line(data=TopSlope,aes(x=date, y=S, colour ="MidSlope"), size=1) +
  xlab("Date") + ylab("S (mm)") + ggtitle("TMWB Model")

## Plot Curve Number Storage at various Ia fractions
p2 <- ggplot() +
  geom_line(data=BotSlopeCN,aes(x=date, y=S, colour ="BotSlope"), size=1) +
  geom_line(data=TopSlopeCN,aes(x=date, y=S, colour ="TopSlope"), size=1) +
  geom_line(data=MidSlopeCN,aes(x=date, y=S, colour ="MidSlope"), size=1) +
  xlab("Date") + ylab("Storage (mm)") + ggtitle("CN Model, IaFrac=0.05")

# rerun lines 105-124 with IaFrac=0.10

p3 <- ggplot() +
  geom_line(data=BotSlopeCN,aes(x=date, y=S, colour ="BotSlope"), size=1) +
  geom_line(data=TopSlopeCN,aes(x=date, y=S, colour ="TopSlope"), size=1) +
  geom_line(data=MidSlopeCN,aes(x=date, y=S, colour ="MidSlope"), size=1) +
  xlab("Date") + ylab("Storage (mm)") + ggtitle("CN Model, IaFrac=0.10")

# rerun lines 105-124 with IaFrac=0.20

p4 <- ggplot() +
  geom_line(data=BotSlopeCN,aes(x=date, y=S, colour ="BotSlope"), size=1) +
  geom_line(data=TopSlopeCN,aes(x=date, y=S, colour ="TopSlope"), size=1) +
  geom_line(data=MidSlopeCN,aes(x=date, y=S, colour ="MidSlope"), size=1) +
  xlab("Date") + ylab("Storage (mm)") + ggtitle("CN Model, IaFrac=0.20")

# rerun lines 105-124 with IaFrac=0.50

p5 <- ggplot() +
  geom_line(data=BotSlopeCN,aes(x=date, y=S, colour ="BotSlope"), size=1) +
  geom_line(data=TopSlopeCN,aes(x=date, y=S, colour ="TopSlope"), size=1) +
  geom_line(data=MidSlopeCN,aes(x=date, y=S, colour ="MidSlope"), size=1) +
  xlab("Date") + ylab("Storage (mm)") + ggtitle("CN Model, IaFrac=0.50")

p1 + p2 + p3 + p4 + p5 + plot_layout(ncol = 3, widths = c(1, 1))

## Homework 2 - calculating QPred

# Running to calculate Qpred for TMWB data
BotSlope$Qpred=NA
BotSlope$Qpred[1]=0
attach(BotSlope)
#fcres=.3
for (t in 2:length(date)){
  S[t]=S[t-1]+Excess[t]     
  Qpred[t]=fcres*S[t]
  S[t]=S[t]-Qpred[t]
}
BotSlope$S=S
BotSlope$Qpred=Qpred 
detach(BotSlope) 
rm(list=c("Qpred","S"))

NSE_BotSlope <- NSE(BotSlope$Qmm,BotSlope$Qpred)

p6 <- ggplot() +
  geom_line(data=BotSlope,aes(x=date, y=Qpred, colour ="BotSlope"), size=1) +
  xlab("Date") + ylab("Predicted Discharge (mm)") + ggtitle("Predicted Discharge using TMWB")

# Running to calculate Qpred for CN data
BotSlopeCN$Qpred=NA
BotSlopeCN$Qpred[1]=0
attach(BotSlopeCN)
#fcres=.3
for (t in 2:length(date)){
  S[t]=S[t-1]+Excess[t]     
  Qpred[t]=fcres*S[t]
  S[t]=S[t]-Qpred[t]
}
BotSlopeCN$S=S
BotSlopeCN$Qpred=Qpred 
detach(BotSlopeCN) 
rm(list=c("Qpred","S"))

NSE_BotSlopeCN <- NSE(BotSlopeCN$Qmm,BotSlopeCN$Qpred)

p7 <- ggplot() +
  geom_line(data=BotSlopeCN,aes(x=date, y=Qpred, colour ="BotSlopeCN"), size=1) +
  xlab("Date") + ylab("Predicted Discharge (mm)") + ggtitle("Predicted Discharge using CN")

p6 + p7 + plot_layout(ncol = 1, widths = c(1, 1))


# Homework 3

source("https://raw.githubusercontent.com/SELoomis/BSE4304/main/Lab05_NewYork")
# This will calculate everything for Rochester

# Plot predicted vs observed data for two watersheds
p8 <- ggplot() +
  geom_line(data=TMWB2new,aes(x=date, y=Qpred, colour ="Predicted"), size=1) +
  geom_line(data=TMWB2new,aes(x=date, y=Qmm, colour ="Observed"), size=1) +
  xlab("Date") + ylab("Discharge (mm)") + ggtitle("Rochester: Predicted vs Observed Discharge")

p9 <- ggplot() +
  geom_line(data=TMWBnew,aes(x=date, y=Qmm, colour ="Observed"), size=1) +
  geom_line(data=TMWBnew,aes(x=date, y=Qpred, colour ="Predicted"), size=1) +
  xlab("Date") + ylab("Discharge (mm)") + ggtitle("Lick Run: Predicted vs Observed Discharge")

p8 + p9 + plot_layout(ncol = 1, widths = c(1, 1))

NSE_LickRun <- NSE(TMWBnew$Qmm,TMWBnew$Qpred)
NSE_Rochester <- NSE(TMWB2new$Qmm,TMWB2new$Qpred)

# Calculate Lick Run's snowmelt
attach(TMWBnew)
SNO_Energy=SnowMelt(date, P, MaxTemp-3, MinTemp-3, myflowgage$declat, 
                       slope = atan(0/100),
                       aspect = 0*(pi/180), tempHt = 1, windHt = 2, groundAlbedo = 0.25,
                       SurfEmissiv = 0.95, windSp = 2, forest = 0, startingSnowDepth_m = 0,
                       startingSnowDensity_kg_m3=450)
detach(TMWBnew)

# Plot snow water equivalent for two watersheds

ggplot() +
  geom_line(data=SNO_Energy_ny,aes(x=Date, y=SnowWaterEq_mm, colour ="Rochester"), size=1) +
  geom_line(data=SNO_Energy,aes(x=Date, y=SnowWaterEq_mm, colour ="Lick Run"), size=1) +
  xlab("Date") + ylab("SWE (mm)") + ggtitle("Snow Water Equivalent for two watersheds")

# Plot snowmelt for two watersheds

ggplot() +
  geom_line(data=SNO_Energy_ny,aes(x=Date, y=SnowMelt_mm, colour ="Rochester"), size=1) +
  geom_line(data=SNO_Energy,aes(x=Date, y=SnowMelt_mm, colour ="Lick Run"), size=1) +
  xlab("Date") + ylab("SWE (mm)") + ggtitle("Snowmelt for two watersheds")
