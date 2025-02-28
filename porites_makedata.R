rm(list=ls())
library(terra)

#source("~/Dropbox/keys_rdaforest_december2024/FL_Keys_rasters_dec2024/smooth.krig.R")
source("~/Dropbox/keys_rdaforest_december2024/FL_Keys_rasters_dec2024/lonlat2raster.R")

ll=load("../porites/Porites_env_rasters_noSmooth_fixVar_gebco.RData")
Depth=env$depth

ll=load("../agaricia/Agaricia_post_jan2025.RData")

r1=lonlat2raster(XY,rasters.post)
el=read.table("../Porites_pops.txt",header=T)
latlon=el[,c("Latitude","Longitude")]
names(latlon)=c("lat","lon")

env=terra::extract(r1,latlon[,2:1])[,-1]
env$depth=Depth

save(XY,latlon,rasters.post,env,file="Porites_post_jan2025.RData")

