library(rerddap)
library(dplyr)
library(ggplot2)
library(viridis)
library(terra)
library(raster)
#setwd("~/Dropbox/FL_Keys_rasters_dec2024")
ll=load("../FL_Keys_rasters_dec2024/SERC_krigs_sst_kd490_prepost.RData")

datasets=ed_search_adv(, query="Chl",minTime = "1997-01-01T00:00:00Z",
                       maxTime="2019-12-31T00:00:00Z", maxLat = 26, minLat = 24, minLon = -83.7, maxLon = -79.5)
print(datasets$info[,1],n=90)
datasets$info[9,]
dataset_id=as.character(datasets$info[9,2])
# erdMH1chlamday
data <- griddap(dataset_id,
                longitude = c(-83.7, -79.5), 
                latitude = c(24, 26), 
                time = c("2003-01-16T00:00:00Z", "2007-12-31"), 
                fields = "all"
)

data.tab=data$data
data.tab$year=sub("-.+","",data.tab$time)
head(data.tab)
save(data.tab,file="chla_old.RData")

# --- fixing stupid rounding problem, resulting in uneven grid
# (physically rewriting lon, lat with constant step)

data.tab$longitude=round(data.tab$longitude,4)
data.tab$latitude=round(data.tab$latitude,4)
#table(is.na(data.tab$k490))

lons=unique(data.tab$longitude)
dd=c()
for(i in 2:length(lons)){
  dd=c(dd,lons[i]-lons[i-1])
}
table(dd)
Step=round(mean(as.numeric(names(table(dd)))),4)
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
head(dll)

Ks=dll %>% 
  group_by(lon,lat,year) %>%
  summarize(chl_median=median(chlorophyll,na.rm=T),chl_hi=quantile(chlorophyll,0.9,na.rm=T),chl_lo=quantile(chlorophyll,0.1,na.rm=T))
Ks$chl_range=Ks$chl_hi-Ks$chl_lo
head(Ks)
Ks1=Ks %>%
  group_by(lon,lat) %>%
  summarize(chl_median=mean(chl_median,na.rm=T),chl_range=mean(chl_range,na.rm=T))

ggplot(Ks1,aes(lon,lat,color=chl_range))+geom_point(pch=15,cex=1)+
  coord_equal()+theme_bw()+scale_color_viridis()

Ks.r=rast(Ks1)
plot(Ks.r)

sst_grid=terra::extract(x = Ks.r, y = XY, method="bilinear" )[-1]

rasters.pre=cbind(rasters.pre,sst_grid)
plot(rasters.pre$chl_median~rasters.pre$k490_median)
#save(latlon.b.pre,rasters.pre,file="latlon_rasters_sst_b_pre.RData")

rr=cbind(XY,sst_grid)
ggplot(rr,aes(lon,lat,color=Chl_range))+geom_point(pch=15,cex=1)+
  coord_equal()+theme_bw()+scale_color_viridis()

#------------ post k490
data <- griddap(dataset_id,
                longitude = c(-83.7, -79.5), 
                latitude = c(24, 26), 
                time = c("2013-01-01", "2022-05-16T00:00:00Z"), 
                fields = "all"
)
data.tab=data$data
data.tab$year=sub("-.+","",data.tab$time)
head(data.tab)
save(data.tab,file="Chl_new.RData")

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
  summarize(Chl_med=median(chlorophyll,na.rm=T),Chl_hi=quantile(chlorophyll,0.9,na.rm=T),Chl_lo=quantile(chlorophyll,0.1,na.rm=T))

ggplot(Ks,aes(lon,lat,color=Chl_med))+geom_point(pch=15,cex=1)+
  coord_equal()+theme_bw()+scale_color_viridis()

Ks$Chl_range=Ks$Chl_hi-Ks$Chl_lo
names(Ks)

Ks.r=rast(Ks)
plot(Ks.r)

sst_grid=terra::extract(x = Ks.r, y = XY, method="bilinear" )[-1]

rasters.post=cbind(rasters.post,sst_grid)
#save(latlon.b.pre,rasters.pre,file="latlon_rasters_sst_b_pre.RData")

rr=cbind(XY,sst_grid)
ggplot(rr,aes(lon,lat,color=Chl_med))+geom_point(pch=15,cex=1)+
  coord_equal()+theme_bw()+scale_color_viridis()
save(rasters.pre,rasters.post,XY,file="SERC_krigs_sst_kd490_Chl.RData")

