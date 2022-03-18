### BSE 4304 - Lab 07
### Sarah Loomis
### 18 March 2022 

objects()
rm(list=objects())
#
# Build a working directory for this weeks lab and change working dir
# Note you might have to specify the path explicitly  as some 
# computers in the lab were not working correctly, to do this go to
# >Session>Set Working Directory
dir.create("~/Week07")
setwd("~/Week07")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(EcoHydRology,curl,httr,rnoaa,raster,shapefiles,rgdal,elevatr, soilDB, ggplot2, patchwork)

# Get our gold standard flow data from USGS 0205551460 `
myflowgage_id="01415460"
myflowgage=get_usgs_gage(myflowgage_id,begin_date = "2015-01-01",
                         end_date = "2019-01-01")

# We want Q in mm/day for the basin
myflowgage$flowdata$Qmm = myflowgage$flowdata$flow/myflowgage$area/10^3

# convert lat/long to UTM projection
trunc((180+myflowgage$declon)/6+1) # gives you the utm zone
proj4_utm = paste0("+proj=utm +zone=", trunc((180+myflowgage$declon)/6+1), " +datum=WGS84 +units=m +no_defs")
print(proj4_utm)

# Lat/Lon (_ll) is much easier!
proj4_ll = "+proj=longlat"

# Now we will build our proj4strings which define our “Coordinate 
# Reference Systems” or CRS in future geographic manipulations. 
crs_ll=CRS(proj4_ll)
crs_utm=CRS(proj4_utm)
print(crs_ll)
print(crs_utm)
# see figure in lab - zone 17 is correct

myflowgage$area   # area in km2

# If the watershed was square, which it is not, the size would be the 
# square root of the area. Also, the gage/pour point is not in the center
# so we will search around the gage.

# Build sp point for USGS gage location, in both _ll and _utm
latlon <- cbind(myflowgage$declon,myflowgage$declat)
myflowgage$gagepoint_ll <- SpatialPoints(latlon)
# note there's no crs -run next line
proj4string(myflowgage$gagepoint_ll)=proj4_ll
myflowgage$gagepoint_utm=spTransform(myflowgage$gagepoint_ll,crs_utm)
# 4.1 mil m north of equator, 594k m north of western boundary of utm zone

# Open up maps.google.com to guesstimate area/lengths
url=paste0("https://www.google.com/maps/@",
           myflowgage$declat,",",myflowgage$declon,",18z")
browseURL(url)

# We are going to over-estimate our area
sqrt(myflowgage$area)   # guestimating square watershed
# For our search we are going to multiply the area by 8 and
# to get the distance
sqrt(myflowgage$area*8)

searchlength=sqrt(myflowgage$area*8)*1000 
pourpoint=SpatialPoints(myflowgage$gagepoint_utm@coords,proj4string = crs_utm)

# add points to this box
bboxpts=myflowgage$gagepoint_utm@coords # just the pour point
bboxpts=rbind(bboxpts,bboxpts+searchlength)
bboxpts=rbind(bboxpts,bboxpts-searchlength)
bboxpts=rbind(bboxpts,c(min(bboxpts[,1]),max(bboxpts[,2])))
bboxpts=rbind(bboxpts,c(max(bboxpts[,1]),min(bboxpts[,2])))
bboxpts=SpatialPoints(bboxpts,proj4string = crs_utm)
# From Lab04, get your DEM
mydem=get_aws_terrain(locations=bboxpts@coords, 
                      z = 12, prj = proj4_utm,src ="aws",expand=1)
#z=12 talks about how high above the terrain you're looking at

res(mydem) # the resolution
plot(mydem)
plot(bboxpts,add=T)
plot(pourpoint,add=T,col="red")

# run terminal stuff here

old_path <- Sys.getenv("PATH")
old_path
Sys.setenv(PATH = paste(old_path,
                        paste0(Sys.getenv("HOME"),"/TauDEM/bin"), 
                        sep = ":"))
system("mpirun aread8")

writeRaster(mydem,filename = "mydem.tif",overwrite=T)

##### all this copied from https://hydrology.usu.edu/taudem/taudem5/TauDEMRScript.txt and edited #####
z=raster("mydem.tif")
plot(z)

# Pitremove
system("mpiexec -n 8 pitremove -z mydem.tif -fel mydemfel.tif")
fel=raster("mydemfel.tif")
plot(fel)


# D8 flow directions
system("mpiexec -n 8 d8flowdir -p mydemp.tif -sd8 mydemsd8.tif -fel mydemfel.tif",show.output.on.console=F,invisible=F)
p=raster("mydemp.tif")
plot(p)
sd8=raster("mydemsd8.tif")
plot(sd8)

# Contributing area
system("mpiexec -n 8 aread8 -p mydemp.tif -ad8 mydemad8.tif")
ad8=raster("mydemad8.tif")
plot(log(ad8))
#zoom(log(ad8))


# Grid Network 
system("mpiexec -n 8 gridnet -p mydemp.tif -gord mydemgord.tif -plen mydemplen.tif -tlen mydemtlen.tif")
gord=raster("mydemgord.tif")
plot(gord)
#zoom(gord)

# DInf flow directions
system("mpiexec -n 8 dinfflowdir -ang mydemang.tif -slp mydemslp.tif -fel mydemfel.tif",show.output.on.console=F,invisible=F)
ang=raster("mydemang.tif")
plot(ang)
slp=raster("mydemslp.tif")
plot(slp)


# Dinf contributing area
system("mpiexec -n 8 areadinf -ang mydemang.tif -sca mydemsca.tif")
sca=raster("mydemsca.tif")
plot(log(sca))
#zoom(log(sca))

# Threshold
system("mpiexec -n 8 threshold -ssa mydemad8.tif -src mydemsrc.tif -thresh 2000")
src=raster("mydemsrc.tif")
plot(src)
plot(pourpoint, add=T)
zoom(src)
plot(pourpoint, add=T)
zoom(src)
plot(pourpoint,add=T)
# until you get the exact pixel


# a quick R function to write a shapefile
makeshape.r=function(sname="shape",n=1)
{
  xy=locator(n=n)
  points(xy)
  
  #Point
  dd <- data.frame(Id=1:n,X=xy$x,Y=xy$y)
  ddTable <- data.frame(Id=c(1),Name=paste("Outlet",1:n,sep=""))
  ddShapefile <- convert.to.shapefile(dd, ddTable, "Id", 1)
  write.shapefile(ddShapefile, sname, arcgis=T)
}

makeshape.r("ApproxOutlets")
# asks for pour point

# Move Outlets
system("mpiexec -n 8 moveoutletstostrm -p mydemp.tif -src mydemsrc.tif -o ApproxOutlets.shp -om Outlet.shp")
outpt=read.shp("Outlet.shp")
approxpt=read.shp("ApproxOutlets.shp")

plot(src)
points(outpt$shp[2],outpt$shp[3],pch=19,col=2)
points(approxpt$shp[2],approxpt$shp[3],pch=19,col=4)

#zoom(src)
# moves point to the center of the raster cell


# Contributing area upstream of outlet
system("mpiexec -n 8 aread8 -p mydemp.tif -o Outlet.shp -ad8 mydemssa.tif")
ssa=raster("mydemssa.tif")
plot(ssa) 


# Threshold
system("mpiexec -n 8 threshold -ssa mydemssa.tif -src mydemsrc1.tif -thresh 4000")
src1=raster("mydemsrc1.tif")
plot(src1)
#zoom(src1)

# Stream Reach and Watershed
system("mpiexec -n 8 streamnet -fel mydemfel.tif -p mydemp.tif -ad8 mydemad8.tif -src mydemsrc1.tif -o Outlet.shp -ord mydemord.tif -tree mydemtree.txt -coord mydemcoord.txt -net mydemnet.shp -w mydemw.tif")
plot(raster("mydemord.tif"),add=T)
#zoom(raster("mydemord.tif"))
plot(raster("mydemw.tif"),add=T)

# Plot streams using stream order as width
snet=read.shapefile("mydemnet")
ns=length(snet$shp$shp)
for(i in 1:ns)
{
  lines(snet$shp$shp[[i]]$points,lwd=snet$dbf$dbf$Order[i])
}

# Peuker Douglas stream definition
system("mpiexec -n 8 peukerdouglas -fel mydemfel.tif -ss mydemss.tif")
ss=raster("mydemss.tif")
plot(ss, add=T)
#zoom(ss)

#  Accumulating candidate stream source cells
system("mpiexec -n 8 aread8 -p mydemp.tif -o Outlet.shp -ad8 mydemssa.tif -wg mydemss.tif")
ssa=raster("mydemssa.tif")
plot(ssa)

#  Drop Analysis
system("mpiexec -n 8 dropanalysis -p mydemp.tif -fel mydemfel.tif -ad8 mydemad8.tif -ssa mydemssa.tif -drp mydemdrp.txt -o Outlet.shp -par 5 500 10 0")

# Stream raster by threshold
system("mpiexec -n 8 threshold -ssa mydemssa.tif -src mydemsrc2.tif -thresh 4000")
plot(raster("mydemsrc2.tif"))

# Stream network
system("mpiexec -n 8 streamnet -fel mydemfel.tif -p mydemp.tif -ad8 mydemad8.tif -src mydemsrc2.tif -ord mydemord2.tif -tree mydemtree2.dat -coord mydemcoord2.dat -net mydemnet2.shp -w mydemw2.tif -o Outlet.shp",show.output.on.console=F,invisible=F)

plot(raster("mydemw2.tif"))
snet=read.shapefile("mydemnet2")
ns=length(snet$shp$shp)
for(i in 1:ns)
{
  lines(snet$shp$shp[[i]]$points,lwd=snet$dbf$dbf$Order[i])
}
#zoom(raster("mydemw2.tif"))

# Wetness Index
system("mpiexec -n 8 slopearearatio -slp mydemslp.tif -sca mydemsca.tif -sar mydemsar.tif", show.output.on.console=F, invisible=F)
sar=raster("mydemsar.tif")
wi=sar
wi[,]=-log(sar[,])
plot(wi)

# Distance Down
system("mpiexec -n 8 dinfdistdown -ang mydemang.tif -fel mydemfel.tif -src mydemsrc2.tif -m ave v -dd mydemdd.tif",show.output.on.console=F,invisible=F)
plot(raster("mydemdd.tif"))


##### Masking ######

mydemw=raster("mydemw.tif")
mybasinmask=trim(mydemw,padding=2)
mybasindem=crop(mydem,mybasinmask)
mybasindem=mask(mybasindem,mybasinmask)
plot(mybasindem)

# Wetness Index
mybasinslp=crop(slp,mybasinmask)
mybasinslp=mask(mybasinslp,mybasinmask)
plot(mybasinslp)

mybasinsca=crop(sca,mybasinmask)
mybasinsca=mask(mybasinsca,mybasinmask)
plot(mybasinsca)

TI = log( (mybasinsca+1)/(mybasinslp+0.00001) )
plot(TI)
#zoom(TI)

pacman::p_load(classInt)
nTIclass=5 #number of TI classes, currently equal area, can adjust method various ways e.g., classIntervals(v, n = nTIclass, style = "jenks")
v=values(TI)
v=v[!is.na(v)]
brks.qt = classIntervals(v, n = nTIclass, style = "quantile")$brks #length nTIclass+1 of just the numeric breakpoints
TIC = cut(TI, breaks=brks.qt, include.lowest = T, right=T)
#
# A series of plots to show all of the components
#

par(mfrow = c(2, 2))
plot(TI, xlab="Long. Distance (m)", ylab="Lat. Distance (m)", main="T.I.")
plot(TIC, xlab="Long. Distance (m)", ylab="Lat. Distance (m)", main="T.I. Classes")
plot(mybasinsca, xlab="Long. Distance (m)", ylab="Lat. Distance (m)", main="FilledDEM - DEM")
plot(mybasinslp, xlab="Long. Distance (m)", ylab="Lat. Distance (m)", main="DEM")
dev.off()

## Homework 3 ##

writeOGR(mydemw_poly,dsn=".",layer="mydemw",driver="ESRI Shapefile", overwrite_layer=TRUE)
# We will use this ESRI shape file, a zipped group of it, to download 
# our soil extent from the WebSoilSurvey Website
zip("mydemw.zip",list.files(pattern="mydemw[:.:]"))
# Download to your local machine mydemw.zip from the "Files" tab

# Open the WebSoilSurvey site to: 
browseURL("https://websoilsurvey.sc.egov.usda.gov/App/WebSoilSurvey.aspx")
# "Creat AOI from a zipped shapefile"
# Open "Download Soils Data" Tab
# "Create Download Link" in lower right hand corner
# Right-Click on download link and "Copy Link Address" and 
# paste into a url object
url="https://websoilsurvey.sc.egov.usda.gov/DSD/Download/AOI/ioaxosgpdblyri42104ptypw/wss_aoi_2022-03-17_18-47-16.zip"
download.file(url,"wss_aoi_2022-03-17_18-47-16.zip") # change to match
unzip("wss_aoi_2022-03-17_18-47-16.zip") # change to match

mysoil=readOGR("wss_aoi_2022-03-17_18-47-16/spatial/soilmu_a_aoi.shp")    
ggplot() + 
  geom_polygon(data = mysoil, aes(x = long, y = lat, group = group), colour = "black", fill = NA)

url="https://www.mrlc.gov/downloads/sciweb1/shared/mrlc/download-tool/NLCD_znIE3Q9eLuxpbPIeLWUY.zip"
download.file(url,"NLCD_znIE3Q9eLuxpbPIeLWUY.zip") # change to match
unzip("NLCD_znIE3Q9eLuxpbPIeLWUY.zip") # change to match

mynlcd=raster("NLCD_2019_Land_Cover_L48_20210604_znIE3Q9eLuxpbPIeLWUY.tiff")
#mybasinnlcd=crop(mynlcd,mybasinmask)
plot(mynlcd)

par(mfrow = c(1, 2))
plot(TIC, xlab="Long. Distance (m)", ylab="Lat. Distance (m)", main="T.I. Classes")
plot(mybasinslp, xlab="Long. Distance (m)", ylab="Lat. Distance (m)", main="DEM")
dev.off()

