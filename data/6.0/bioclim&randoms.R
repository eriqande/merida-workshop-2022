library(raster)
library(sp)
library(rgdal)

setwd("~/Documents/CurrentProjects/UCLA_UCSC/BGP/WorldClim/wc2.0_30s_bio")
list.files()
bio09=raster("wc2.0_bio_30s_09.tif")
bio15=raster("wc2.0_bio_30s_15.tif")
bio01=raster("wc2.0_bio_30s_01.tif")
bio13=raster("wc2.0_bio_30s_13.tif")
bio11=raster("wc2.0_bio_30s_11.tif")
bio18=raster("wc2.0_bio_30s_18.tif")
srtm=raster("../srtm/srtm.tif")

##Lets look at one raster
# Plot raster object
plot(srtm,main="World climate raster of BIO09, Precipitation of warmest quarter")

##Combine rasters into one stack
srtmstack=stack(srtm)
bioclim=stack(bio09,bio15,bio01,bio13,bio11,bio18)

##extract for specific lat/long
coord<-meta %>% select(Lat,Long)
srtm_pabu<-extract(srtmstack,coord,fun=mean)
bioclim_pabu<-extract(bioclim,coord,fun=mean)
##Create random sites for Willow Flycatcher breeding range
##Obviously change the path to wherever your shape file is
pabu=readOGR(dsn = path.expand("~/Dropbox/BGP/genoscape_maps/shapefiles/PABU/"), layer = "Passerina_ciris")
summary(pabu)
breed = subset(pabu,SEASONAL=2)
proj4string(pabu)=CRS("+proj=longlat +datum=WGS84")
proj4string(breed)=CRS("+proj=longlat +datum=WGS84")
plot(breed)
##Limit the Climate variables to JUST the breeding range of Willow Flycatcher
bioclimw=crop(bioclim,pabu)
bioclimpabu=mask(bioclimw,pabu)
summary(bioclimpabu)
plot(bioclimpabu)
srtmstackw=crop(srtmstack,pabu)
srtmstackpabu=mask(srtmstackw,pabu)
summary(srtmstackpabu)

##Did it work? Let's extract a single raster from our stack 
bio09pabu <- raster(bioclimpabu, layer=1)
srtmpabu <- raster(srtmstackpabu, layer=1)
all<-stack(bioclimw,srtmstackw)
all<-stack(bioclimw,srtmstackw)

##AND plot bio 18 now that everywhere but the breeding range is masked
plot(bio09pabu,
     main="BIO09 limited to Willow Flycatcher breeding range")
plot(srtmpabu,
     main="SRTM limited to Willow Flycatcher breeding range")

#Sample random sites only within breeding range, rename the columns and then print out your random sites
randomsitesP=sampleRandom(bioclimpabu,size=100000,na.rm=TRUE,xy=TRUE,sp=TRUE)
randomsitesPS=sampleRandom(srtmpabu,size=100000,na.rm=TRUE,xy=TRUE,sp=TRUE)
?raster
head(randomsitesP)
head(randomsitesPS)
randomsitesP=as.data.frame(randomsitesP)
colnames(randomsitesP)=c("X","Y","BIO9","BIO15","BIO1","BIO13")
write.table(randomsitesP,file="pabugridCB.txt",row.names=FALSE,quote=F,sep='\t')

randomsitesPS=as.data.frame(randomsitesPS)
colnames(randomsitesPS)=c("X","Y","SRTM")
write.table(randomsitesPS,file="pabugrid_srtmCB.txt",row.names=FALSE,quote=F,sep='\t')