library(rerddap)
library(dplyr)
library(ggplot2)
library(viridis)
library(terra)
library(raster)
setwd("~/Dropbox/FL_Keys_rasters_dec2024")
ll=load("SERC_krigs_sst.RData")

datasets=ed_search_adv(, query="kd490",minTime = "1995-01-01T00:00:00Z",
                       maxTime="2020-02-01T00:00:00Z", maxLat = 26, minLat = 24, minLon = -83.7, maxLon = -79.5)

dataset_id=as.character(datasets$info[1,2])
# erdMH1kd4901day
data <- griddap(dataset_id,
                longitude = c(-83.7, -79.5), 
                latitude = c(24, 26), 
                time = c("2008-05-01T00:00:00Z", "2019-05-01T00:00:00Z"), 
                fields = "all"
)
data.tab=data$data
data.tab$year=sub("-.+","",data.tab$time)
head(data.tab)
save(data.tab,file="k490_may2008_may2019.RData")

# --- fixing stupid rounding problem, resulting in uneven grid
# (physically rewriting lon, lat with constant step)

data.tab$longitude=round(data.tab$longitude,5)
data.tab$latitude=round(data.tab$latitude,5)
#table(is.na(data.tab$k490))

lons=unique(data.tab$longitude)
dd=c()
for(i in 2:length(lons)){
  dd=c(dd,lons[i]-lons[i-1])
}
table(dd)
Step=round(mean(as.numeric(names(table(dd)))),5)
lats=unique(data.tab$latitude)

length(lats)
# 49
# odd number, starting with median
latmed=round(median(lats),5)
latup=seq(latmed+Step,length.out=floor(length(lats)/2),by=Step)
latdn=rev(seq(latmed-Step,length.out=floor(length(lats)/2),by=-Step))
relat=rev(c(latdn,latmed,latup))
length(relat)
range(relat)
plot(relat~lats)

length(lons)
# 101
# odd number, starting with median
lonmed=round(median(lons),5)
lonup=seq(lonmed+Step,length.out=floor(length(lons)/2),by=Step)
londn=rev(seq(lonmed-Step,length.out=floor(length(lons)/2),by=-Step))
relon=c(londn,lonmed,lonup)
length(relon)
range(relon)
plot(relon~lons)

lontab=data.frame(longitude=lons,lon=relon)
lattab=data.frame(latitude=lats,lat=relat)

dl=merge(data.tab,lontab,by="longitude")
dll=merge(dl,lattab,by="latitude")
dll=dll[,-c(1:2)]

Ks=dll %>% 
  group_by(lon,lat) %>%
  summarize(k490_med=median(k490,na.rm=T),k490_hi=quantile(k490,0.9,na.rm=T),k490_lo=quantile(k490,0.1,na.rm=T))

ggplot(Ks,aes(lon,lat,color=k490_med))+geom_point(pch=15,cex=1)+
  coord_equal()+theme_bw()+scale_color_viridis()

Ks$k490_range=Ks$k490_hi-Ks$k490_lo
names(Ks)

Ks.r=rast(Ks)
plot(Ks.r)

sst_grid=terra::extract(x = Ks.r, y = XY, method="bilinear" )[-1]

rasters.pre=cbind(rasters.pre,sst_grid)
#save(latlon.b.pre,rasters.pre,file="latlon_rasters_sst_b_pre.RData")

rr=cbind(XY,sst_grid)
ggplot(rr,aes(lon,lat,color=k490_range))+geom_point(pch=15,cex=1)+
  coord_equal()+theme_bw()+scale_color_viridis()

#------------ post k490
data <- griddap(dataset_id,
                longitude = c(-83.7, -79.5), 
                latitude = c(24, 26), 
                time = c("2008-01-01", "2019-12-31"), 
                fields = "all"
)
data.tab=data$data
data.tab$year=sub("-.+","",data.tab$time)
head(data.tab)
save(data.tab,file="k490_new.RData")

# --- fixing stupid rounding problem, resulting in uneven grid
# (physically rewriting lon, lat with constant step)

data.tab$longitude=round(data.tab$longitude,5)
data.tab$latitude=round(data.tab$latitude,5)
#table(is.na(data.tab$k490))

lons=unique(data.tab$longitude)
dd=c()
for(i in 2:length(lons)){
  dd=c(dd,lons[i]-lons[i-1])
}
table(dd)
Step=round(mean(as.numeric(names(table(dd)))),5)
lats=unique(data.tab$latitude)

length(lats)
# 49
# odd number, starting with median
latmed=round(median(lats),5)
latup=seq(latmed+Step,length.out=floor(length(lats)/2),by=Step)
latdn=rev(seq(latmed-Step,length.out=floor(length(lats)/2),by=-Step))
relat=rev(c(latdn,latmed,latup))
length(relat)
range(relat)
plot(relat~lats)

length(lons)
# 101
# odd number, starting with median
lonmed=round(median(lons),5)
lonup=seq(lonmed+Step,length.out=floor(length(lons)/2),by=Step)
londn=rev(seq(lonmed-Step,length.out=floor(length(lons)/2),by=-Step))
relon=c(londn,lonmed,lonup)
length(relon)
range(relon)
plot(relon~lons)

lontab=data.frame(longitude=lons,lon=relon)
lattab=data.frame(latitude=lats,lat=relat)

dl=merge(data.tab,lontab,by="longitude")
dll=merge(dl,lattab,by="latitude")
dll=dll[,-c(1:2)]

Ks=dll %>% 
  group_by(lon,lat) %>%
  summarize(k490_med=median(k490,na.rm=T),k490_hi=quantile(k490,0.9,na.rm=T),k490_lo=quantile(k490,0.1,na.rm=T))

ggplot(Ks,aes(lon,lat,color=k490_med))+geom_point(pch=15,cex=1)+
  coord_equal()+theme_bw()+scale_color_viridis()

Ks$k490_range=Ks$k490_hi-Ks$k490_lo
names(Ks)

Ks.r=rast(Ks)
plot(Ks.r)

sst_grid=terra::extract(x = Ks.r, y = XY, method="bilinear" )[-1]

rasters.post=cbind(rasters.post,sst_grid)
#save(latlon.b.pre,rasters.pre,file="latlon_rasters_sst_b_pre.RData")

rr=cbind(XY,sst_grid)
ggplot(rr,aes(lon,lat,color=k490_med))+geom_point(pch=15,cex=1)+
  coord_equal()+theme_bw()+scale_color_viridis()

save(rasters.pre,rasters.post,XY,file="SERC_krigs_sst_kd490.RData")

