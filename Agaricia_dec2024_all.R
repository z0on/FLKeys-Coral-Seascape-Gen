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
dim(XY)
dim(rasters.post)
hist(distinct(env)$dhw_max)

plot(lonlat2raster(XY,data.frame(rasters.post[,c("SRP.B_median","SRP.B_range","DIN.TP_range")])),bty="n")

# admixture proportions between 3 clusters
admix=read.csv("Agaricia_admix.csv")[,-1]

# genetic distances
IBS=as.matrix(read.table("Agaricia_nozoox.ibsMat"))
dim(IBS)
samples=row.names(env)

# adding sample names to everything
dimnames(IBS)=list(samples,samples)
row.names(admix)=row.names(latlon)=row.names(env)=samples

# ------------ choosing subpop to work on
 pop="all"
# 
#  #choices0=row.names(env)[which(admix[,pop]>=0.95)]
#  #which(samples[choices] %in% c("aa212","aa527"))
#  choices=scan(paste("Ag",pop,"_nozoox60",sep=""),what="c")
#  choices=paste("aa",sub(".bam|_.+","",choices),sep="")
#  # reading lineage-specific IBS, making sure its samples match to main IBS
#  IBS.p=as.matrix(read.table(paste("Ag",pop,"_nozoox60.ibsMat",sep="")))
#  dimnames(IBS.p)=list(choices,choices)
#  latlon.p=latlon[choices,]
#  admix.p=admix[choices,]
#  env.p=env[choices,]
#  admix.cov=admix.p[,1:3]
#  admix.cov=admix.cov[,!(names(admix.cov) %in% pop)]

IBS.p=IBS
latlon.p=latlon
env.p=env
admix.p=admix
choices=row.names(env.p)

mm=plot_nice_map(latlon.p[choices,],mapdata=mapdata,margin=40000,cols="grey50",pch=3,size=2)

# converting to UTM coordinates
epsg.code=epsg.maker(mean(range(latlon.p$lon)),mean(range(latlon.p$lat)))
latlon.p.utm=latlon2UTM(latlon.p, epsg.code)
row.names(latlon.p.utm)=row.names(IBS.p)

# IBD plot
plot(as.dist(IBS.p)~dist(latlon.p.utm),pch=19,cex=0.5,col=rgb(0,0,0,0.1))

protest(capscale(dist(latlon.p.utm)~1),capscale(as.dist(IBS.p)~1))
# all: R=0.14, p=0.02

plot(capscale(IBS.p~1+Condition(as.matrix(cbind(latlon.p.utm)))),scaling=1)

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
#admix.cov=admix.cov[goods,]
env.p=env.p[goods,]

# ---- modification to ordination command when analyzing individual lineages,
#      to DISREGARD admixture with other lineages as a force shaping genetic distances
#      (conditional PCoA)

# computing conditional PCoA

ord.all=capscale(as.dist(IBS.p)~1+Condition(as.matrix(cbind(latlon.p.utm))))
plot(ord.all,scaling=1)
plot(ord.all$CA$eig/sum(ord.all$CA$eig))

# GF model including spatials
gf=makeGF(ord.all,cbind(env.p),ntrees=1500,pcs2keep=c(1:20))
# which MDSses are predicted?
gf$result
# with dis2land
# MDS1       MDS4      MDS10      MDS20 
# 0.39838883 0.16012844 0.24812893 0.01420351 

# MDS1       MDS4      MDS10      MDS20      MDS43 
# 0.39710204 0.14899493 0.24509023 0.00184261 0.01891732 
tokeep=11

# plot importances (R2)
imps=data.frame(importance_RDAforest(gf,ord.all))
# total R2
sum(imps)
# 0.332
names(imps)="R2"
imps$var=row.names(imps)
imps$var=factor(imps$var,levels=imps$var[order(imps$R2)])
pdf(file=paste(pop,"_importances_simple.pdf",sep=""),height=6,width=4)
ggplot(imps,aes(var,R2))+geom_bar(stat="identity")+coord_flip()+theme_bw()
dev.off()

pdf(file=paste(pop,"_turnovers_simple.pdf",sep=""),height=6,width=4)
plot_gf_turnovers(gf,imps$var[c(1:4)])
dev.off()


mm=mtrySelJack(Y=IBS.p,X=cbind(env.p),nreps=30,covariates=cbind(latlon.p.utm),top.pcs=tokeep,importance.cutoff = 0.05,prop.positive.cutoff = 0.5)
mm$goodvars
# [1] "depth"
mm=Reselect(mm,0.5,0.05)
#  "depth"      "k490_range"
# dis2land loses
save(mm,gf,file=paste(pop,"_mtrySel.RData",sep=""))
load(paste(pop,"_mtrySel.RData",sep=""))

# boxplot of importance differences at different mtry
pop="all"
pdf(file=paste(pop,"_mtrySel.pdf",sep=""),height=6,width=5)
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
imps=data.frame(cbind(mm$importances))
names(imps)="R2"
imps$var=row.names(imps)
imps$var=factor(imps$var,levels=imps$var[order(imps$R2)])
ggplot(imps,aes(var,R2))+geom_bar(stat="identity")+coord_flip()+theme_bw()
dev.off()

#--------------- jackknifing of mtry-passing variables plus latlon
tokeep=11
sb.old=ordinationJackknife(Y=IBS.p,X=cbind(env.p)[,mm$goodvars],
                       covariates=cbind(latlon.p.utm),
                       newX=cbind(XY,rasters.pre),
                       nreps=20,top.pcs=tokeep,extra=1)
sb=ordinationJackknife(Y=IBS.p,X=env.p[,mm$goodvars],
                       covariates=cbind(latlon.p.utm),
                       newX=cbind(XY,rasters.post),
                       nreps=20,top.pcs=tokeep,extra=0.1)
sum(sb$median.importance)
# 0.326

save(sb,sb.old,mm,gf,file=paste(pop,"_rdaForest01.RData",sep=""))
load(paste(pop,"_rdaForest.RData",sep=""))

OFFS=gen_offset_oj(X=sb.old,Y=sb,sx=XY,sy=XY)
head(OFFS)
dim(OFFS)
# plotting the offset as a raster on the map (square coordinates):
pdf(file=paste(pop,"_offset01.pdf",sep=""),width=6,height=4)
# # simple way
plot(lonlat2raster(OFFS[,1:2],values=data.frame(cbind(offset=OFFS$offset))),col=(map.pal("viridis")))
cols=viridis_colors(OFFS$offset)
plot_nice_map(OFFS[,1:2],mapdata=mapdata,margin=40000,cols=cols,pch=15,size=0.17)
dev.off()

names(env)
env$sst_range=env$sst_max-env$sst_min
pairs(env[,c("sst_max","sst_min","sst_range","dhw_max","TEMP.B_median","TEMP.B_range")])
cor(env[,c("sst_max","sst_min","sst_range","dhw_max","TEMP.B_median","TEMP.B_range")])
cor(env$sst_range,env$sst_)

pdf(file=paste(pop,"_turnovers.pdf",sep=""),height=3.5,width=3)
plot_turnover(sb,rasters.post[sb$goodrows,],mm$goodvars[1])
plot_turnover(sb,rasters.post[sb$goodrows,],mm$goodvars[2])
dev.off()


#------ PLOTTING
# sb_old=sb
# sb=sb_new
# sb=sb_old
pop="all"
load(paste(pop,"_rdaForest.RData",sep=""))
pdf(file=paste(pop,"_importance01.pdf",sep=""),height=1.5,width=3)
important=names(sb$median.importance)
ggplot(sb$all.importances[sb$all.importances$variable %in% important,],aes(variable,importance))+
  geom_boxplot(outlier.shape = NA)+
  coord_flip()+
#  ylim(0,0.2)+
  theme_bw()
# what is the proportion of variation explained, total?
dev.off()


phys=names(sb$median.importance)[grep("k490|DSIG",names(sb$median.importance))]
therm=names(sb$median.importance)[grep("TEMP|sst|dhw",names(sb$median.importance))]
chem=names(sb$median.importance)[!(names(sb$median.importance) %in% c(phys,therm,"depth"))]

Depth=sb$median.importance["depth"]
Phys=sum(sb$median.importance[phys])
Therm=sum(sb$median.importance[therm])
Chem=sum(sb$median.importance[chem])
meds=data.frame(cbind(predictor=c("depth","Therm","Phys","Chem")))
meds$sum.R2=c(Depth,Therm,Phys,Chem)
meds$predictor=factor(meds$predictor,levels=rev(c("depth","Therm","Phys","Chem")))

pdf(file=paste(pop,"_importance_summed.pdf",sep=""),height=2,width=2.5)
ggplot(meds,aes(predictor,sum.R2))+
  geom_bar(stat="identity")+
#  ylim(0,0.22)+
  coord_flip()+
  #  ylim(0,0.2)+
  theme_bw()
dev.off()
sum(meds$sum.R2)
#0.326

#-----------------
pop="all"
load(paste(pop,"_rdaForest01.RData",sep=""))
turnovers=sb$predictions.turnover
raster.vars=colnames(turnovers)

table(sb$goodrows)
xy2=XY[sb$goodrows,]
dirs=sb$predictions.direct



ras2=rasters.post[which(sb$goodrows),]
bests=names(sb$median.importance)[c(1:2)]

pdf(file=paste(pop,"_bestVars.pdf",sep=""),height=6,width=12)
plot(lonlat2raster(xy2,ras2[,bests]),bty="n")
dev.off()


bests=names(sb$median.importance)[c(1,2)]
pdf(file=paste(pop,"_maps_extra01.pdf",sep=""),width=8,height=6)
#plot_adaptation(dirs,ras2[,bests],xy2,color.scheme="100",scal=10,size=15,cex=0.8,jitscale=0.05,rangeExp=2.5)
pa=plot_adaptation(dirs,ras2[,bests],xy2,nclust=20,
                   cluster.guide = turnovers,
                   pcs2plot=c(2,3),
                   cluster.merge.level = 0.4, 
                   cluster.labels=T,ordinate=F,pc.jitter=0.005,
                   color.scheme="011",scal=10,cex=0.5,jitscale=0.025,rangeExp=2,
                   lighten=0.3)
# pa0=plot_adaptation(dirs,ras2[,bests],xy2,nclust=12,
#                    cluster.merge.level = 0.4, 
#                    cluster.labels=T,color.scheme="100",
#                    scal=8,cex=0.4,rangeExp=1.5,
#                    jitscale=0.025,lighten=0.3)
plot_nice_map(xy2,size=0.35,mapdata=mapdata,map.on.top=T,cols=pa$colors,pch=15)
plot_nice_map(xy2,size=0.31,mapdata=mapdata,map.on.top=T,cols=pa$colors,pch=15,overlay.points = latlon[,2:1],pch.points = 3,cols.points ="grey40")
dev.off()

names(env)
