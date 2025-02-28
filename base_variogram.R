library(dplyr)
library(ggplot2)
library(viridis)
library(terra)
library(raster)
library(vegan)
library(gstat)

ll=load("rasters_keys_dec2024_noNA.RData")

names(rasters.pre)

sats=rasters.pre[,c(39,40,41,44,45,48,49)]
# save(sats,file="satellite.RData")
# rows_with_na <- apply(sats, 1, function(x) any(is.na(x)))
# bads=which(rows_with_na)
# xy=XY[-bads,]
# sats=sats[-bads,]
# 
# rasters.pre=rasters.pre[-bads,]
# rasters.post=rasters.pre[-bads,]
# XY=XY[-bads,]
# save(rasters.pre,rasters.post,XY,file="rasters_keys_dec2024_noNA.RData")
# 
# ss=data.frame(cbind(xy,sats))
# ggplot(ss,aes(lon,lat,color=depth))+geom_point()

rr=rda(sats~1,scale=T)
plot(rr,scaling=1,display="sites")
plot(rr$CA$eig/sum(rr$CA$eig))
scc=data.frame(scores(rr,scaling=1,display='sites'))


ggs=data.frame(cbind(XY,scc))
ggplot(ggs,aes(lon,lat,color=PC1))+geom_point()+scale_color_viridis()
sat.pc.scores=scc
save(XY,colors,sat.pc.scores,file="sat.pc.scores.RData")


variogram_data <- data.frame(cbind(XY,z=sat.pc.scores$PC1))
coordinates(variogram_data) <- ~lon + lat
proj4string(variogram_data) <- CRS("+proj=longlat +datum=WGS84")
vdata_utm <- spTransform(variogram_data, CRS("+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"))
dd=extent(vdata_utm)
# max distance from corner to corner
maxd=dist(rbind(c(dd@xmin,dd@ymin),c(dd@xmax,dd@ymax)))
maxd
# 350224

# Calculate variogram and fit model
variogram_empirical <- variogram(z ~ 1, data = vdata_utm,cutoff = 30000, width = 500)
plot(variogram_empirical)
variogram_model <- fit.variogram(variogram_empirical, model = vgm(psill = 1, model = "Mat", range = 10000, nugget = 0.1,kappa=1))
pdf("matern_variogram_satellitePC1.pdf",width=5,height=4)
plot(variogram_empirical, variogram_model, main = "Fitted Variogram")
dev.off()
save(variogram_empirical, variogram_model,file="variogram_model_satPC1.RData")


     