library(rerddap)
library(dplyr)
library(ggplot2)
library(viridis)
library(terra)
library(raster)
source("~/Dropbox/keys_rdaforest_december2024/FL_Keys_rasters_dec2024/lonlat2raster.R")

ll=load("~/Dropbox/keys_rdaforest_december2024/FL_Keys_rasters_dec2024/SERC_krigs_sst_kd490_prepost.RData")
names(rasters.pre)
library(sf)

depth=raster("~/Dropbox/GEBCO_bathymetry/gebco_2024_n26.1154_s24.1128_w-84.6087_e-79.7214.tif")
# depth2=raster("~/Dropbox/noaa_bathymetry.tiff")
#depth=raster("bathymetry.grd")

d_grid=terra::extract(x = depth, y = XY)

rasters.pre$depth=(-1)*d_grid
rasters.post$depth=(-1)*d_grid
rpre=lonlat2raster(XY,rasters.pre)
dim(rpre)
plot(rpre[[28:32]])
rpost=lonlat2raster(XY,rasters.post)

save(rasters.pre,rasters.post,XY,file="rasters_keys_jan2025_gebcoDepth.RData")









ll=load("Agaricia_envdata_clean_oldDepth_dec2024.RData")
ll=load("rasters_keys_dec2024_gebcoDepth.RData")
save(rasters.pre,rasters.post,XY,env,latlon,file="agaricia_clean_dec2024_gebcoDepth.RData")

names(rasters.pre)
sats=rasters.pre[,28:32]
rr=rda(sats~1,scale=T)
# plot(rr,scaling=1,display="sites")
# plot(rr$CA$eig/sum(rr$CA$eig))
scc=data.frame(scores(rr,scaling=1,display='sites'))

# making a colorful map plot based on pcs (skip if just need variogram)
b= scc$PC2
g = scc$PC1-scc$PC2
r = scc$PC1
r <- (r - min(r))/(max(r) - min(r)) * 255
g <- (g - min(g))/(max(g) - min(g)) * 255
b <- (b - min(b))/(max(b) - min(b)) * 255
colors=rgb(r,g,b,max=255)
plot(XY,pch=16,cex=0.3,col=colors,asp=1)
# plotting PC1
ggplot(cbind(XY,PC1=scc$PC1),aes(lon,lat,color=PC1))+
  geom_point(cex=0.03)+scale_color_viridis()+coord_equal()+theme_minimal()


# converting data to UTM coordinates
variogram_data <- data.frame(cbind(XY,z=scc$PC1))
coordinates(variogram_data) <- ~lon + lat
proj4string(variogram_data) <- CRS("+proj=longlat +datum=WGS84")
vdata_utm <- spTransform(variogram_data, CRS("+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"))
dd=extent(vdata_utm)
# max distance from corner to corner
maxd=dist(rbind(c(dd@xmin,dd@ymin),c(dd@xmax,dd@ymax)))
maxd
# 350224

# empirical variogram
variogram_empirical = variogram(z ~ 1, data = vdata_utm,cutoff = 40000, width = 500)
plot(variogram_empirical)

# fitting a Matern model. It is quite robust to vgm() parameter settings.
variogram_model <- fit.variogram(variogram_empirical, model = vgm(psill = 1, model = "Mat", range = 10000, nugget = 10,kappa=1))
plot(variogram_empirical, variogram_model, main = "Fitted Variogram")
save(variogram_empirical, variogram_model,file="variogram_model_satPC1.RData")

# the object variogram_model can be passed to smooth.krig() as argument variogram_model, 
# to force it not to fit its own model but use the supplied one.


