### Lab 11
### Sarah Loomis
### April 15 2022

if (!require("pacman")) install.packages("pacman")
pacman::p_load(httr,EcoHydRology,curl,data.table,multisensi)
#
setwd("~")
objects()
rm(list=objects())
#
dir.create("~/Week11")
setwd("~/Week11/")
list.files(all.files = T)
objects()   # Should be empty.

T <- seq(from = 5, to = 365, by = 5)

##### PET_fromTemp sensitivity analysis #####

#PET_fromTemp(Jday, Tmax_C, Tmin_C, lat_radians, 
#             AvgT = (Tmax_C + Tmin_C)/2, albedo = 0.18, TerrestEmiss = 0.97, 
#             aspect = 0, slope = 0, forest = 0, PTconstant=1.26,
#             AEparams=list(vp=NULL, opt="linear"))

PET_fromTemp2 <- function (Jday, Tmax_C, Tmin_C, lat_radians, AvgT = (Tmax_C + Tmin_C)/2, albedo = 0.18, TerrestEmiss = 0.97, aspect = 0, slope = 0, forest = 0, PTconstant=1.26, AEparams=list(vp=NULL, opt="linear"))
{
  cloudiness <- EstCloudiness(Tmax_C, Tmin_C)
  DailyRad <- NetRad(lat_radians, Jday, Tmax_C, Tmin_C, albedo, forest, slope, aspect, AvgT, cloudiness, TerrestEmiss, AvgT, AEparams=AEparams)
  potentialET <- PTpet(DailyRad, AvgT, PTconstant)
  potentialET[which(potentialET < 0)] <- 0
  potentialET[which(Tmax_C == -999 | Tmin_C == -999)] <- (-999)
  return(potentialET)
}

PET2 <- function(X, t = T) {
  out <- matrix(nrow = nrow(X), ncol = length(t), NA)
  for (i in 1:nrow(X)) {
    out[i, ] <- PET_fromTemp2(Jday=t, 
                                    Tmax_C=X$Tmax_C[i],
                                    Tmin_C=(X$Tmax_C[i]-X$Trange[i]),
                                    lat_radians=X$lat_radians[i],
                                    AvgT = (X$Tmax_C[i] + (X$Tmax_C[i]-X$Trange[i]))/2, 
                                    aspect=X$aspect[i],
                                    slope=X$slope[i])
  }
  out <- as.data.frame(out)
  names(out) <- paste("t", t, sep = "")
  return(out)
}

n <- 10
set.seed(1234)
X <- data.frame(Tmax_C = runif(n, min = 1, max = 50), 
                Trange = runif(n, min = 1, max = 10),
                lat_radians = runif(n,min=0,max=pi/3),
                slope=runif(n, min=0, max=0.2),
                aspect=runif(n, min=0, max=pi*2))
Y <- PET2(X)
#par(cex.axis = 0.7, cex.lab = 0.8)
#plot(T, Y[1, ], type = "l", xlab = "Time", ylab = "Population size",
#     ylim = c(0, max(Y)))
#for (i in 2:n) {
#  lines(T, Y[i, ], type = "l", col = i)
#}

X <- expand.grid(Tmax_C = c(1,30,50), 
                 Trange = c(1,5,10), 
                 lat_radians = c(0, pi/4, pi/3),
                 slope = c(0.1,0.2,0.3),
                 aspect = c(0.1,.5,1.0)) #not a random distribution this time
Y <- PET2(X) ## this part can be performed outside R if necessary
PET.seq <- multisensi(design=X, model=Y, reduction=NULL, center=FALSE)
## [*] Analysis + Sensitivity Indices

library(multisensi)
PET.seq.fast <- multisensi(design = fast99, model = PET2,
                           center = FALSE, reduction = NULL, analysis = analysis.sensitivity,
                           design.args=list( factors=c("Tmax_C","Trange","lat_radians", "slope","aspect"), n=1000, q = "qunif",
                                             q.arg = list(list(min=1, max=50), list(min=1, max=10),
                                                          list(min = 0, max = pi/3), list(min=0,max=0.2), list(min=0,max=pi*2))),
                           analysis.args=list(keep.outputs=FALSE))


PET.pca <- multisensi(design=X, model=Y, reduction=basis.ACP, scale=FALSE)
summary(PET.pca, digits = 2)


plot(PET.pca, graph = 1)
plot(PET.pca, graph = 2)
plot(PET.pca, graph = 3)


##### NetRad sensitivity analysis #####

#NetRad(lat, Jday, Tx, Tn, albedo = 0.18, forest = 0, slope = 0, 
#       aspect = 0, airtemp = (Tn+Tx)/2, cloudiness = "Estimate", 
#       surfemissivity = 0.97, surftemp = (Tn+Tx)/2, units = "kJm2d", 
#       AEparams=list(vp=NULL, opt="linear"))

NetRad2 <- function(X, t = T) {
  out <- matrix(nrow = nrow(X), ncol = length(t), NA)
  for (i in 1:nrow(X)) {
    out[i, ] <- EcoHydRology::NetRad(lat=X$lat[i], 
                                    Tx=X$Tx[i], 
                                    Tn=(X$Tx[i]-X$Trange[i]), 
                                    albedo=X$albedo[i],
                                    slope=X$slope[i],
                                    aspect=X$aspect[i],
                                    Jday=t,  
                                    units="kJm2d")
  }
  out <- as.data.frame(out)
  names(out) <- paste("t", t, sep = "")
  return(out)
}

n <- 10
set.seed(1234)
X <- data.frame(lat = runif(n, min = 0, max = pi/3), 
                Tx = runif(n, min = 1, max = 40), 
                Trange = runif(n, min = 1, max = 10),
                albedo=runif(n, min=0, max=1),
                slope=runif(n, min=0, max=0.2),
                aspect=runif(n, min=0, max=pi*2))
Y <- NetRad2(X)
#par(cex.axis = 0.7, cex.lab = 0.8)
#plot(T, Y[1, ], type = "l", xlab = "Time (days)", ylab = "Daily Net Radiation (kJ/m2d)",
#     ylim = c(0, max(Y)))
#for (i in 2:n) {
#  lines(T, Y[i, ], type = "l", col = i)
#}

library(multisensi)
X <- expand.grid(Tx = c(1,30,50), 
                 Trange = c(1,5,10), 
                 lat = c(0, pi/4, pi/3),
                 slope = c(0.1,0.2,0.3),
                 aspect = c(0.1,.5,1.0),
                 albedo = c(0.1, 0.3, 0.8)) #not a random distribution this time
Y <- NetRad2(X) ## this part can be performed outside R if necessary
NetRad.seq <- multisensi(design=X, model=Y, reduction=NULL, center=FALSE)

NetRad.seq.fast <- multisensi(design = fast99, model = NetRad2,
                             center = FALSE, reduction = NULL, analysis = analysis.sensitivity,
                             design.args=list( factors=c("lat", "Tx","Trange","albedo", "slope", "aspect"), n=1000, q = "qunif",
                                               q.arg = list(list(min=0, max=pi/3), list(min=1, max=40),
                                                            list(min = 1, max = 10), list(min=0, max=1), list(min=0, max=0.2),
                                                            list(min=0, max=pi*2))),
                             analysis.args=list(keep.outputs=FALSE))
## [*] Design
## [*] Response simulation
## [*] Analysis + Sensitivity Indices
NetRad.pca <- multisensi(design=X, model=Y, reduction=basis.ACP, scale=FALSE)
summary(NetRad.pca, digits = 2)

plot(NetRad.pca, graph = 1)
plot(NetRad.pca, graph = 2)
plot(NetRad.pca, graph = 3)
