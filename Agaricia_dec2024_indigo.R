rm(list=ls())
library(RDAforest)
library(rnaturalearth)
library(rnaturalearthdata)
library(terra)
library(viridis)
library(sf)
library(spatstat)

mapdata= ne_countries(scale="large",continent="north america")
source("~/Dropbox/keys_rdaforest_december2024/FL_Keys_rasters_dec2024/smooth.krig.R")
source("~/Dropbox/keys_rdaforest_december2024/FL_Keys_rasters_dec2024/lonlat2raster.R")

ll=load("~/Dropbox/keys_rdaforest_december2024/FL_Keys_rasters_dec2024/rasters_keys_feb2025.RData")
ll=load("~/Dropbox/keys_rdaforest_december2024/FL_Keys_rasters_dec2024/Agaricia_env_feb2025.RData")
viridis_colors <- function(input_vector,...) {
  normalized_values <- scales::rescale(input_vector,to=c(0,1))  # Normalize to [0,1]
  colors <- viridis(100,...)  # Generate 100 Viridis colors for smooth mapping
  color_indices <- round(normalized_values * (length(colors) - 1)) + 1  # Scale to color indices
  return(colors[color_indices])  # Assign colors based on values
}

# admixture proportions between 3 clusters
admix=read.csv("Agaricia_admix.csv")[,-1]
samples=paste("aa",sub(".bam","",admix$Sample),sep="")
# genetic distances
IBS=as.matrix(read.table("Agaricia_nozoox.ibsMat"))
dim(IBS)

tt=grep("sst|TEMP",names(rasters.post))
pdf("thermal_rasters.pdf",width=7,height=4)
plot(lonlat2raster(XY, rasters.post[,tt]))
dev.off()

# genetic distances
IBS=as.matrix(read.table("Agaricia_nozoox.ibsMat"))
dim(IBS)

# adding sample names to everything
dimnames(IBS)=list(samples,samples)
row.names(admix)=row.names(latlon)=row.names(env)=samples

# ------------ choosing subpop to work on
 pop="indigo"

 #choices0=row.names(env)[which(admix[,pop]>=0.95)]
 #which(samples[choices] %in% c("aa212","aa527"))
 choices=scan(paste("Ag",pop,"_nozoox90",sep=""),what="c")
 choices=paste("aa",sub(".bam|_.+","",choices),sep="")
 length(choices)
 # 62
 # reading lineage-specific IBS, making sure its samples match to main IBS
 IBS.p=as.matrix(read.table(paste("Ag",pop,"_nozoox90.ibsMat",sep="")))
 dimnames(IBS.p)=list(choices,choices)
 latlon.p=latlon[choices,]
 admix.p=admix[choices,]
 env.p=env[choices,]
 admix.cov=admix.p[,1:3]
 admix.cov=admix.cov[,!(names(admix.cov) %in% pop)]

# IBS.p=IBS
# latlon.p=latlon
# env.p=env1
# admix.p=admix


mm=plot_nice_map(latlon.p,mapdata=mapdata,margin=40000,cols="gold",pch=3,size=2)

# converting to UTM coordinates
epsg.code=epsg.maker(mean(range(latlon.p$lon)),mean(range(latlon.p$lat)))
latlon.p.utm=latlon2UTM(latlon.p, epsg.code)
row.names(latlon.p.utm)=row.names(IBS.p)

# IBD plot
plot(as.dist(IBS.p)~dist(latlon.p.utm),pch=19,cex=0.5,col=rgb(0,0,0,0.1))

protest(capscale(dist(latlon.p.utm)~1),capscale(as.dist(IBS.p)~1))
# all: R=0.14, p=0.02
# yellow: R=0.80, p<0.001
# blue R=0.54, p<0.001
# indigo: R=0.64, p<0.001

plot(capscale(IBS.p~1+Condition(as.matrix(cbind(latlon.p.utm,admix.cov)))),scaling=1)

# ------ hierarchical clustering tree to look for clones, wrong species collections

hc=hclust(as.dist(IBS.p),"ave")
plot(hc,cex=0.5) # clustering of samples by IBS (great to detect clones or closely related individuals)
# there are a few "low-hanging" groups of closely related individuals
# Note: to make sure that two samples are clones, their distance must be compared to the distance
# between genotyping replicates (the same sample genotyped twice)
abline(h=0.01,col="red") # this seems like a reasonable  "low-hang" threshold for calling related groups
#pheatmap::pheatmap(1-IBS.p)

#--------- retaining a single representative of each highly-related group

cuts=cutree(hc,h=0.01)
goods=c();i=1
for (i in unique(cuts)) {
  goods[i]=names(cuts[cuts==i])[1]
}
length(goods)  # how many samples are left?
# subsetting all data for only the retained samples
IBS.p=IBS.p[goods,goods]
latlon.p.utm=latlon.p.utm[goods,]
admix.cov=admix.cov[goods,]
env.p=env.p[goods,]

# ---- modification to ordination command when analyzing individual lineages,
#      to DISREGARD admixture with other lineages as a force shaping genetic distances
#      (conditional PCoA)

# computing conditional PCoA

ord.all=capscale(as.dist(IBS.p)~1+Condition(as.matrix(cbind(latlon.p.utm,admix.cov))))
plot(ord.all,scaling=1)
plot(ord.all$CA$eig/sum(ord.all$CA$eig))

# GF model including spatials
gf=makeGF(ord.all,cbind(env.p,latlon.p.utm,admix.cov),ntrees=1500,pcs2keep=c(1:20))
# which MDSses are predicted?
gf$result
# MDS1       MDS2       MDS3       MDS4       MDS5       MDS6      MDS10 
# 0.72275848 0.37044856 0.13702509 0.20514414 0.20339088 0.18870940 0.05485649 
tokeep=11

# plot importances (R2)
imps=data.frame(importance_RDAforest(gf,ord.all))
# total R2
sum(imps)
# 0.06
names(imps)="R2"
imps$var=row.names(imps)
imps$var=factor(imps$var,levels=imps$var[order(imps$R2)])
pdf(file=paste(pop,"_importances_simple_dl.pdf",sep=""),height=6,width=4)
ggplot(imps,aes(var,R2))+geom_bar(stat="identity")+coord_flip()+theme_bw()
dev.off()

pdf(file=paste(pop,"_turnovers_simple_post.pdf",sep=""),height=6,width=4)
plot_gf_turnovers(gf,imps$var[c(1:4)])
dev.off()


mm=mtrySelJack(Y=IBS.p,X=cbind(env.p),nreps=30,covariates=cbind(latlon.p.utm,admix.cov),top.pcs=tokeep,importance.cutoff = 0.1,prop.positive.cutoff = 0.5)
mm$goodvars
# with dis2land
# [1] "SRP.B_median" "depth"        "k490_median"  "NH4.B_range"  "SRP.B_range"  "DIN.B_median"
# [7] "k490_range"   "DSIGT_range"  "dhw_max"      "dis2land"     "N.P_range"   
># [1] "SRP.B_median" "depth"        "k490_median"  "NH4.B_range"  "SRP.B_range"  "NO2.B_median" "DSIGT_range" 
# [8] "dhw_max" 
save(mm,gf,tokeep,file=paste(pop,"_mtrySel.RData",sep=""))
#load(paste(pop,"_mtrySel.RData",sep=""))

# boxplot of importance differences at different mtry
pdf(file=paste(pop,"_mtrySel_dl.pdf",sep=""),height=6,width=5)
# ggplot(mm$delta,aes(var,values))+
#   geom_boxplot()+
#   coord_flip()+
#   geom_hline(yintercept=0,col="red")

# bar chart of proportion of positive change in response to higher mtry
# good variables would be the ones above the red line
ggplot(mm$prop.positive,aes(var,prop.positive))+
  geom_bar(stat="identity")+
  coord_flip()+
  geom_hline(yintercept=0.5,col="red")
hist(mm$prop.positive$prop.positive)
dev.off()

#--------------- jackknifing of mtry-passing variables plus latlon

sb.old=ordinationJackknife(Y=IBS.p,X=cbind(env.p)[,mm$goodvars],
                       covariates=cbind(latlon.p.utm,admix.cov),
                       newX=cbind(XY,rasters.pre),
                       nreps=25,top.pcs=tokeep,extra=1)
sb=ordinationJackknife(Y=IBS.p,X=cbind(env.p)[,mm$goodvars],
                       covariates=cbind(latlon.p.utm,admix.cov),
                       newX=cbind(XY,rasters.post),
                       nreps=25,top.pcs=tokeep,extra=0.1)


sum(sb.old$median.importance)
# 0.059
sum(sb$median.importance)
# 0.0615

#save(sb.old,sb.new,mm,tokeep,gf,file=paste(pop,"_rdaForest_post.RData",sep=""))
save(sb,sb.old,mm,tokeep,gf,file=paste(pop,"_rdaForest_dl.RData",sep=""))

OFFS=gen_offset_oj(X=sb.old,Y=sb,sx=XY,sy=XY)
head(OFFS)
dim(OFFS)
# plotting the offset as a raster on the map (square coordinates):
pdf(file=paste(pop,"_offset01.pdf",sep=""),width=6,height=4)
# # simple way
plot(lonlat2raster(OFFS[,1:2],values=data.frame(cbind(offset=OFFS$offset))),col=(map.pal("viridis")))
plot_nice_map(OFFS[,1:2],mapdata=mapdata,margin=40000,cols=viridis_colors(OFFS$offset),pch=15,size=0.17)
dev.off()


pdf(file=paste(pop,"_turnovers.pdf",sep=""),height=3.5,width=3)
plot_turnover(sb,rasters.post[sb$goodrows,],"dhw_max")
plot_turnover(sb,rasters.post[sb$goodrows,],mm$goodvars[10])
plot_turnover(sb,rasters.post[sb$goodrows,],mm$goodvars[2])
plot_turnover(sb,rasters.post[sb$goodrows,],mm$goodvars[3])
plot_turnover(sb,rasters.post[sb$goodrows,],mm$goodvars[4])
dev.off()

#------ PLOTTING
# sb=sb_old
pdf(file=paste(pop,"_importance_dl.pdf",sep=""),height=2.5,width=3)
important=names(sb$median.importance)
ggplot(sb$all.importances[sb$all.importances$variable %in% important,],aes(variable,importance))+
  geom_boxplot(outlier.shape = NA)+
  coord_flip()+
#  ylim(0,0.2)+
  theme_bw()
# what is the proportion of variation explained, total?
dev.off()
sum(sb$median.importance)
# I: 0.0615

pop="indigo"
load(paste(pop,"_rdaForest_dl.RData",sep=""))
phys=names(sb$median.importance)[grep("k490|DSIG",names(sb$median.importance))]
therm=names(sb$median.importance)[grep("TEMP|sst|dhw",names(sb$median.importance))]
chem=names(sb$median.importance)[!(names(sb$median.importance) %in% c(phys,therm,"depth"))]

Dis2land=sb$median.importance["dis2land"]
Depth=sb$median.importance["depth"]
Phys=sum(sb$median.importance[phys])
Therm=sum(sb$median.importance[therm])
Chem=sum(sb$median.importance[chem])
meds=data.frame(cbind(predictor=c("dis2land","depth","Therm","Phys","Chem")))
meds$sum.R2=c(Dis2land,Depth,Therm,Phys,Chem)
meds$predictor=factor(meds$predictor,levels=rev(c("depth","dis2land","Therm","Phys","Chem")))

pdf(file=paste(pop,"_importance_summed.pdf",sep=""),height=2,width=2.5)
ggplot(meds,aes(predictor,sum.R2))+
  geom_bar(stat="identity")+
  #  ylim(0,0.22)+
  coord_flip()+
  #  ylim(0,0.2)+
  theme_bw()
dev.off()
sum(meds$sum.R2)
# 0.06

#-----------------
pop="indigo"
ll=load(paste(pop,"_rdaForest01.RData",sep=""))
#ll=load("Agaricia_post_jan2025.RData")

turnovers=sb$predictions.turnover
raster.vars=colnames(turnovers)

table(sb$goodrows)
xy2=XY[sb$goodrows,]
dirs=sb$predictions.direct

ras2=rasters.post[which(sb$goodrows),]
bests=names(sb$median.importance)[c(1:4)]

pdf(file=paste(pop,"_bestVars01.pdf",sep=""),height=6,width=12)
plot(lonlat2raster(xy2,ras2[,bests]),bty="n")
dev.off()

pdf(file=paste(pop,"_maps_dl.pdf",sep=""),width=8,height=6)
#plot_adaptation(dirs,ras2[,bests],xy2,color.scheme="000",scal=8,cex=0.8,jitscale=0.05,rangeExp=1.5)
pa=plot_adaptation(dirs,ras2[,bests],xy2,nclust=20,
                   cluster.guide = turnovers,
                   cluster.merge.level = 0.57, 
                   cluster.labels=T,
                   color.scheme="000",scal=10,cex=0.4,jitscale=0.05,rangeExp=2.5,
                   lighten=0.3)
plot_nice_map(xy2,size=0.31,mapdata=mapdata,map.on.top=T,cols=pa$colors,pch=15)
plot_nice_map(xy2,size=0.31,mapdata=mapdata,map.on.top=T,cols=pa$colors,pch=15,overlay.points = latlon[choices,2:1],pch.points = 3,cols.points ="grey40")
dev.off()

# tortugas: -82.8554,24.6652
# looe key: -81.4059,24.5485
# top marquesas: -82.14913461392415,24.61030167326765
lc=c(-81.4059,24.5485)

reef=c( -82.14913461392415,24.61030167326765)
#plot(lonlat2raster(xy2,sb$predictions.direct))
rf.preds=terra::extract(lonlat2raster(xy2,sb$predictions.direct),data.frame(rbind(reef)),method="simple")[,-1]
sc=adapt_scale(sb$predictions.direct)[2]
# computing distances between future-needed and present-day adaptation
agf=env_mismatch(X=rf.preds,Y=sb,sy=XY,sc=sc)

drf=data.frame(rbind(reef))
names(drf)=c("lon","lat")
ggplot(agf,aes(lon,lat,color=agf$env.mismatch))+
  geom_point(cex=0.3,pch=15)+
  coord_equal()+
  theme_minimal()+
  scale_color_viridis(direction=-1)+
  geom_point(data=data.frame(rbind(reef)),aes(X1,X2),pch=3,size=2,col="red")

# Function to generate Viridis colors for an input vector
cols <- viridis_colors(agf$env.mismatch,direction=-1)
drf=rbind(drf,lc)
plot_nice_map(agf[,1:2],mapdata=mapdata,overlay.points=drf,pch.points=3,size.points = 2,cols.points="red",margin=40000,cols=cols,pch=15,size=0.17)

plot_nice_map(agf[,1:2],mapdata=mapdata,margin=40000,cols=cols,pch=15,size=0.27)


