dir.create("~/Week03Lab03/")
setwd("~/Week03Lab03/")
# The solution to last weeks lab solution: 
# https://drive.google.com/open?id=14fFB1WDGhWdlADl7xLNIbttdZt0-wwTB
# a handy trick in google docs is the export=download parameter option 
# https://docs.google.com/a/vt.edu/uc?id=14fFB1WDGhWdlADl7xLNIbttdZt0-wwTB&export=download
# becomes: 
url="https://docs.google.com/a/vt.edu/uc?id=14fFB1WDGhWdlADl7xLNIbttdZt0-wwTB&export=download"
# This will grab the solution for last weeks Lab02 Homework
download.file(url,"Lab02_HW1_Solution.R")
file.edit("Lab02_HW1_Solution.R")



 SFTmp = 1  # referred to as SFTMP in SWAT input (Table 1)
 bmlt6 = 4.5   # referred to as SMFMX in SWAT input (Table 1)
 bmlt12 = 0.0  # referred to as SMFMN in SWAT input adjusted for season
 Tmlt = SFTmp  # Assumed to be same as SnowFall Temperature
 Tlag = 1  # referred to as TIMP in SWAT input (Table 1)
 TMWB$AvgTemp=(TMWB$MaxTemp + TMWB$MinTemp)/2
   TMWB$bmlt = (bmlt6 + bmlt12)/2 + (bmlt6 - bmlt12)/2 *  sin(2*pi/365*(julian(TMWB$date,origin = as.Date("2000-01-01"))-81))
# Initialize SNO, Tsno as well as the first values of each
 TMWB$SNO = 0  # Snow Depth (mm)
 TMWB$Tsno = 0  # Snow Temp (C)
 TMWB$SNOmlt = 0  # Snow Melt (mm)
 attach(TMWB)
 for (t in 2:length(date)){
  Tsno[t]= Tsno[t-1] * (1.0-Tlag) +  AvgTemp[t] * Tlag
  if(AvgTemp[t] < SFTmp){
    SNO[t]= SNO[t-1] + P[t]
  }  else {
    SNOmlt[t]= bmlt[t] * SNO[t-1] * ((Tsno[t]+MaxTemp[t])/2 - Tmlt) 
    SNOmlt[t]= min(SNOmlt[t],SNO[t-1])
    SNO[t]= SNO[t-1] -SNOmlt[t]
  }
  print(t)
}
 plot(date,SNO,type="l")
 detach(TMWB)
 TMWB$Tsno=Tsno
 TMWB$SNO=SNO
 TMWB$SNOmlt=SNOmlt
 rm(list=c("SNO", "SNOmlt", "Tsno"))

  bmlt12 = 0.0
  bmlt6 = 3
  SFTmp = 7
  TMWB$bmlt = (bmlt6 + bmlt12)/2 + (bmlt6 - bmlt12)/2 *  sin(2*pi/365*(julian(TMWB$date,origin = as.Date("2000-01-01"))-81))
 # Initialize SNO, Tsno as well as the first values of each
  TMWB$SNO = 0  # Snow Depth (mm)
  TMWB$Tsno = 0  # Snow Temp (C)
  TMWB$SNOmlt = 0  # Snow Melt (mm)
  attach(TMWB)
  for (t in 2:length(date)){
   Tsno[t]= Tsno[t-1] * (1.0-Tlag) +  AvgTemp[t] * Tlag
   if(AvgTemp[t] < SFTmp){
     SNO[t]= SNO[t-1] + P[t]
   }  else {
     SNOmlt[t]= bmlt[t] * SNO[t-1] * ((Tsno[t]+MaxTemp[t])/2 - Tmlt) 
     SNOmlt[t]= min(SNOmlt[t],SNO[t-1])
     SNO[t]= SNO[t-1] -SNOmlt[t]
   }
   print(t)
 }
  lines(date,SNO,col="red")
  detach(TMWB)
  TMWB$Tsno=Tsno
  TMWB$SNO=SNO
  TMWB$SNOmlt=SNOmlt
  rm(list=c("SNO", "SNOmlt", "Tsno"))
  
  # notice that there is an Energy Balance based Snow Accumulation 
  # and Melt model in the EcoHydRology package.
   ?SnowMelt
   attach(TMWB)
   SNO_Energy=SnowMelt(date, P, MaxTemp-3, MinTemp-3, myflowgage$declat, 
                        slope = 0,
                        aspect = 0, tempHt = 1, windHt = 2, groundAlbedo = 0.25,
                        SurfEmissiv = 0.95, windSp = 2, forest = 0, startingSnowDepth_m = 0,
                        startingSnowDensity_kg_m3=450)
  # How do we know what units slope and aspect are in? (view function)
   lines(date,SNO_Energy$SnowWaterEq_mm)
   detach(TMWB)
  
TMWB$Albedo=.23
TMWB$Albedo[TMWB$SNO>0]=.95 # if there is snow, albedo is higher
?PET_fromTemp # units in meters
attach(TMWB)
PET=PET_fromTemp(Jday=(1+as.POSIXlt(date)$yday),Tmax_C = MaxTemp,Tmin_C = MinTemp,albedo=Albedo,lat_radians = myflowgage$declat*pi/180) * 1000
# multiply by 1000 for mm
TMWB$PET=PET
plot(date,PET)
detach(TMWB)
rm(list=c("PET"))

myflowgage$FldCap=.45
myflowgage$WiltPt=.15
myflowgage$Z=1000
TMWB$AWC=(myflowgage$FldCap-myflowgage$WiltPt)*myflowgage$Z # 
TMWB$dP = 0 # Initializing Net Precipitation
TMWB$ET = 0 # Initializing ET
TMWB$AW = 0 # Initializing AW
TMWB$Excess = 0 # Initializing Excess


# Loop to calculate AW and Excess
attach(TMWB)
for (t in 2:length(AW)){
 #This is where Net Precipitation is now calculated
  # Do you remember what Net Precip is? Refer to week 2 notes
  ET[t] = min (AW[t-1],PET[t])
  ET[t] = (AW[t-1]/AWC[t-1])*PET[t] # New Model - how much water in soil profile
  if(AvgTemp[t] >= SFTmp){
    dP[t] = P[t] - ET[t] + SNOmlt[t] 
  }  else {
    dP[t] = ET[t]
  }

  # From here onward, everything is the same as Week2’s lab
  if (dP[t]<=0) {
    values<-soildrying(AW[t-1],dP[t],AWC[t])
  } else if((dP[t]>0) & (AW[t-1]+dP[t])<=AWC[t]) {
    values<-soilwetting(AW[t-1],dP[t],AWC[t])
  } else {
    values<-soil_wetting_above_capacity(AW[t-1],dP[t],AWC[t])
  }
  AW[t]<-values[1]
  Excess[t]<-values[2]
  print(t)
}
TMWB$AW=AW
TMWB$Excess=Excess
TMWB$dP=dP
rm(list=c("AW","dP","ET", "Excess"))
detach(TMWB) # IMPORTANT TO DETACH

TMWB$Qpred=NA
TMWB$Qpred[1]=0
TMWB$S=NA
TMWB$S[1]=0
attach(TMWB)
fcres=.3
for (t in 2:length(date)){
  S[t]=S[t-1]+Excess[t]     
  Qpred[t]=fcres*S[t]
  S[t]=S[t]-Qpred[t]
}
TMWB$S=S
TMWB$Qpred=Qpred # UPDATE vector BEFORE DETACHING

#Make a plot that has Qmm, P,and Qpred over time
plot(date,P,col="black")
lines(date,Qmm,type = "l",col="black")
lines(date,Qpred,col="blue")
detach(TMWB) # IMPORTANT TO DETACH
rm(list=c("Qpred","S"))

   
