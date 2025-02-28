library(rerddap)
library(dplyr)
library(ggplot2)
library(viridis)
library(terra)
library(raster)
ll=load("SERC_krigs_bottom_prepost_5k.RData")

datasets1=ed_search_adv(, query='Bleaching',minTime = "2000-01-01T00:00:00Z",
                       maxTime="2010-02-01T00:00:00Z", maxLat = 26, minLat = 24, minLon = -83.7, maxLon = -79.5)
datasets1$title
# datasets2=ed_search_adv(, query="Bleaching",minTime = "2010-01-01T00:00:00Z",
#                        maxTime="2020-02-01T00:00:00Z", maxLat = 26, minLat = 24, minLon = -83.7, maxLon = -79.5)
# 
# datasets2$info[1,]
# 

# Old temperature, 5km resolution
dataset_id <- "NOAA_DHW_monthly"
data <- griddap(dataset_id,
                longitude = c(-83.7, -79.5), 
                latitude = c(24, 26), 
                time = c("2003-01-01", "2007-12-31"), 
                fields = "all"
)
data.tab=data$data[data$data$mask==0,]
data.tab$mask=NULL
data.tab$year=sub("-.+","",data.tab$time)
head(data.tab)

SSTs=data.tab %>% 
  group_by(longitude,latitude,year) %>%
  summarize(sst_min=min(sea_surface_temperature),sst_max=max(sea_surface_temperature),dhw_max=max(sea_surface_temperature_anomaly))

SST=SSTs %>%
  group_by(longitude,latitude) %>%
  summarize(sst_min=mean(sst_min),sst_max=mean(sst_max),dhw_max=mean(dhw_max))

names(SST)[1:2]=c("lon","lat")
SST$lon=round(SST$lon,3)
SST$lat=round(SST$lat,3)

sst.r=rast(SST)
plot(sst.r)

sst_grid=terra::extract(x = sst.r, y = XY, method="bilinear" )[-1]

rasters.pre=cbind(rasters.pre,sst_grid)
#save(latlon.b.pre,rasters.pre,file="latlon_rasters_sst_b_pre.RData")

rr=cbind(XY,sst_grid)
ggplot(rr,aes(lon,lat,color=sst_max))+geom_point(pch=15,cex=1)+coord_equal()
ggplot(rr,aes(lon,lat,color=dhw_max))+geom_point(pch=15,cex=1)+coord_equal()
ggplot(rr,aes(lon,lat,color=sst_min))+geom_point(pch=15,cex=1)+coord_equal()

#------------ post temperature

data <- griddap(dataset_id,
                longitude = c(-83.7, -79.5), 
                latitude = c(24, 26), 
                time = c("2014-01-01", "2019-12-31"), 
                fields = "all"
)
data.tab=data$data[data$data$mask==0,]
data.tab$mask=NULL
data.tab$year=sub("-.+","",data.tab$time)
head(data.tab)

SSTs=data.tab %>% 
  group_by(longitude,latitude,year) %>%
  summarize(sst_min=min(sea_surface_temperature),sst_max=max(sea_surface_temperature),dhw_max=max(sea_surface_temperature_anomaly))

SST=SSTs %>%
  group_by(longitude,latitude) %>%
  summarize(sst_min=mean(sst_min),sst_max=mean(sst_max),dhw_max=mean(dhw_max))

names(SST)[1:2]=c("lon","lat")
SST$lon=round(SST$lon,3)
SST$lat=round(SST$lat,3)

sst.r=rast(SST)
plot(sst.r)

sst_grid=terra::extract(x = sst.r, y = XY, method="bilinear" )[-1]

rasters.post=cbind(rasters.post,sst_grid)

rr=cbind(XY,sst_grid)
ggplot(rr,aes(lon,lat,color=sst_max))+geom_point(pch=15,cex=1)+coord_equal()
ggplot(rr,aes(lon,lat,color=dhw_max))+geom_point(pch=15,cex=1)+coord_equal()
ggplot(rr,aes(lon,lat,color=sst_min))+geom_point(pch=15,cex=1)+coord_equal()


save(rasters.pre,rasters.post,XY,file="SERC_krigs_sst_prepost.RData")

