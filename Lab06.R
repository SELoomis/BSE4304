### BSE 4304 - Lab 06
### Sarah Loomis
### 5 March 2022 

rm(list=objects())
setwd("~/Lab06")
laburl="https://raw.githubusercontent.com/vtdrfuka/BSE5304_2022/main/Lab05/lab05sol"
download.file(laburl,"Lab05Sol.R")
file.edit("Lab05Sol.R")
# Run this script to get all the base data
if (!require("pacman")) install.packages("pacman")
pacman::p_load(elevatr,soilDB,rgdal,raster, patchwork)
pacman::p_load(EcoHydRology,rnoaa,curl,httr,ggplot2)

##### Homework 1: DEoptim for TMWB and CN models #####
CNmodeldf=modeldata

## For CN model
f <- function (x) {
  CNopt=x[1]
  IaOpt=x[2]
  func_DAWCopt=x[3]
  func_zopt=x[4]
  fnc_fcresopt=x[5]
  CNmodelnew=CNmodel(CNmodeldf =CNmodeldf,CNavg = CNopt,IaFrac = IaOpt,func_DAWC = func_DAWCopt,
                     func_z = func_zopt, fnc_fcres = fnc_fcresopt)
  NSE_CN=1-NSE(CNmodelnew$Qmm,CNmodelnew$Qpred) 
  return(NSE_CN)
}

lower <- c(35,.01,0.1,500,0.1)
upper <- c(99,.25,0.3,1500,0.5)

set.seed(1234)
DEoptim(f, lower, upper,control = DEoptim.control(itermax=40))
# optimized at CNopt=95.761657, IaOpt=0.019359, func_DAWCopt=0.277739, func_zopt=1487.310340, and fnc_fcresopt=0.340022
# NSE_CN = 1-0.434062 = 0.565938

CNmodelnew=CNmodel(CNmodeldf =CNmodeldf,CNavg = 95.761657,IaFrac = 0.019359,func_DAWC = 0.277739,
                   func_z = 1487.310340, fnc_fcres = 0.340022)

## For TMWB model
g <- function (x) {
  fcresopt=x[1]
  SFTmpopt=x[2]
  Tlagopt=x[3]
  AWCvalopt=x[4]
  TMWBnew=TMWBmodel(TMWB =TMWB, fcres = fcresopt, SFTmp=SFTmpopt, Tlag=Tlagopt, AWCval=AWCvalopt)
  1-NSE(TMWBnew$Qmm,TMWBnew$Qpred)  
}

lower <- c(0.1,-5,0,150)
upper <- c(0.5,5,1,350)

set.seed(1234)
outDEoptim <- DEoptim(g, lower, upper,control = DEoptim.control(itermax=40))
# optimized at fcresopt=0.318802, SFTmpopt=4.487581, Tlagopt=0.777663, and AWCvalopt=150.016889
# NSE_TMWB = 1-0.695141 = 0.304859

TMWBnew=TMWBmodel(TMWB =TMWB, fcres = 0.318802, SFTmp=4.487581, Tlag=0.777663, AWCval=150.016889)

p1 <- ggplot() +
  geom_line(data=TMWBnew,aes(x=date, y=Qpred), size=1) +
  xlab("Date") + ylab("Discharge (mm)") + ggtitle("TMWB model")
p2 <- ggplot() +
  geom_line(data=CNmodelnew,aes(x=date, y=Qpred), size=1) +
  xlab("Date") + ylab("Discharge (mm)") + ggtitle("CN model")
p1 + p2 + plot_layout(ncol = 1, widths = c(1, 1))


##### Homework 2: VSA-CN model #####

# Initialize the TI Class objects from top to bottom of slope
pacman::p_load(lubridate, data.table)
BasinCN_JO=CNmodelnew[(month(CNmodelnew$date) > 5 
                      & month(CNmodelnew$date) < 11),]
attach(BasinCN_JO)

h <- function (x) {
  Sest=x
  NSE(Qmm,dP^2/(dP+Sest))
}
optimize(h, c(50,500), tol = 0.0001,maximum = TRUE)$maximum
detach(BasinCN_JO)
Sest=153.8569 # this value is based on the result of the optimization

nTIclass=5
VSAsol=data.table(WetClass=seq(from=nTIclass,to=1),
                  As=seq(1:nTIclass)*(1/nTIclass),Wetfrac=(1/nTIclass))
VSAsol[,sSratio:=2*(sqrt(1-shift(As))-sqrt(1-As))/Wetfrac-1]
VSAsol$sSratio[1]=2*(sqrt(1-0)-sqrt(1-VSAsol$As[1]))/VSAsol$Wetfrac[1]-1
VSAsol[,sigma:=Sest*sSratio]
VSAsol[,CN:=25400/(sigma+254)]

TIC05=modeldata
TIC05 = CNmodel(CNmodeldf = TIC05, CNavg=VSAsol$CN[1],
                func_DAWC=.3,IaFrac=0.05,
                func_z=1000,fnc_fcres=.3)

TIC04=modeldata
TIC04$P=TIC05$Excess+TIC04$P
TIC04 = CNmodel(CNmodeldf = TIC04, CNavg=VSAsol$CN[2],
                func_DAWC=.3,IaFrac=0.05,
                func_z=1000,fnc_fcres=.3)

TIC03=modeldata
TIC03$P=TIC04$Excess+TIC03$P
TIC03 = CNmodel(CNmodeldf = TIC03, CNavg=VSAsol$CN[3],
                func_DAWC=.3,IaFrac=0.05,
                func_z=1000,fnc_fcres=.3)

TIC02=modeldata
TIC02$P=TIC03$Excess+TIC02$P
TIC02 = CNmodel(CNmodeldf = TIC02, CNavg=VSAsol$CN[4],
                func_DAWC=.3,IaFrac=0.05,
                func_z=1000,fnc_fcres=.3)

TIC01=modeldata
TIC01$P=TIC02$Excess+TIC01$P
TIC01 = CNmodel(CNmodeldf = TIC01, CNavg=VSAsol$CN[5],
                func_DAWC=.3,IaFrac=0.05,
                func_z=1000,fnc_fcres=.3)

TIC01avg=mean(TIC01$Excess)
TIC02avg=mean(TIC02$Excess)
TIC03avg=mean(TIC03$Excess)
TIC04avg=mean(TIC04$Excess)
TIC05avg=mean(TIC05$Excess)

ggplot() +
  geom_line(data=TIC01,aes(x=date, y=Qpred, colour ="TIC01"), size=1) +
  geom_line(data=TIC02,aes(x=date, y=Qpred, colour ="TIC02"), size=1) +
  geom_line(data=TIC03,aes(x=date, y=Qpred, colour ="TIC03"), size=1) +
  geom_line(data=TIC04,aes(x=date, y=Qpred, colour ="TIC04"), size=1) +
  geom_line(data=TIC05,aes(x=date, y=Qpred, colour ="TIC05"), size=1) +
  xlab("Date") + ylab("Discharge (mm)") + ggtitle("VSA discharge along hillslope")


##### Homework 3: Change in AW #####
TIC01_AWmin=min(TIC01$AW)
TIC01_AWavg=mean(TIC01$AW)
TIC01_AWmax=max(TIC01$AW)

TIC02_AWmin=min(TIC02$AW)
TIC02_AWavg=mean(TIC02$AW)
TIC02_AWmax=max(TIC02$AW)

TIC03_AWmin=min(TIC03$AW)
TIC03_AWavg=mean(TIC03$AW)
TIC03_AWmax=max(TIC03$AW)

TIC04_AWmin=min(TIC04$AW)
TIC04_AWavg=mean(TIC04$AW)
TIC04_AWmax=max(TIC04$AW)

TIC05_AWmin=min(TIC05$AW)
TIC05_AWavg=mean(TIC05$AW)
TIC05_AWmax=max(TIC05$AW)

ggplot() +
  geom_line(data=TIC01,aes(x=date, y=AW, colour ="TIC01"), size=1) +
  geom_line(data=TIC02,aes(x=date, y=AW, colour ="TIC02"), size=1) +
  geom_line(data=TIC03,aes(x=date, y=AW, colour ="TIC03"), size=1) +
  geom_line(data=TIC04,aes(x=date, y=AW, colour ="TIC04"), size=1) +
  geom_line(data=TIC05,aes(x=date, y=AW, colour ="TIC05"), size=1) +
  xlab("Date") + ylab("Available Water (mm)") + ggtitle("AW among TI classes")
