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


# genetic distances
IBS=as.matrix(read.table("Agaricia_nozoox.ibsMat"))
dim(IBS)

# adding sample names to everything
dimnames(IBS)=list(samples,samples)
row.names(admix)=row.names(latlon)=row.names(env)=samples

# ------------ choosing subpop to work on
 pop="blue"

 #choices0=row.names(env)[which(admix[,pop]>=0.95)]
 #which(samples[choices] %in% c("aa212","aa527"))
 choices=scan(paste("Ag",pop,"_nozoox60",sep=""),what="c")
 choices=paste("aa",sub(".bam|_.+","",choices),sep="")
 length(choices)
 # 49
 # reading lineage-specific IBS, making sure its samples match to main IBS
 IBS.p=as.matrix(read.table(paste("Ag",pop,"_nozoox60.ibsMat",sep="")))
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


mm=plot_nice_map(latlon.p[choices,],mapdata=mapdata,margin=40000,cols="skyblue",pch=3,size=2)

# converting to UTM coordinates
epsg.code=epsg.maker(mean(range(latlon.p$lon)),mean(range(latlon.p$lat)))
latlon.p.utm=latlon2UTM(latlon.p, epsg.code)
row.names(latlon.p.utm)=row.names(IBS.p)

# IBD plot
plot(as.dist(IBS.p)~dist(latlon.p.utm),pch=19,cex=0.5,col=rgb(0,0,0,0.1))

protest(capscale(dist(latlon.p.utm)~1),capscale(as.dist(IBS.p)~1))
# blue: R=0.54, p<0.001

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
gf=makeGF(ord.all,cbind(env.p),ntrees=1500,pcs2keep=c(1:20))
# which MDSses are predicted?
gf$result
# with dis2land
# MDS1       MDS2       MDS3       MDS4       MDS5      MDS11 
# 0.57223174 0.52246002 0.64997935 0.33799426 0.07610098 0.05004582
# only env
# MDS1       MDS2       MDS3       MDS4       MDS5      MDS11 
# 0.58171923 0.51678204 0.65323521 0.32091784 0.11796295 0.04407585 
# with admix and latlon (admix is major)
# MDS1       MDS2       MDS3       MDS4       MDS5 
# 0.64435972 0.54335857 0.63688812 0.36207437 0.06607663
# keeping first 6 MDSes
tokeep=12

# plot importances (R2)
imps=data.frame(importance_RDAforest(gf,ord.all))
# total R2
sum(imps)
# 0.171 with admix 
# 0.16 without
names(imps)="R2"
imps$var=row.names(imps)
imps$var=factor(imps$var,levels=imps$var[order(imps$R2)])
pdf(file=paste(pop,"_importances_simple_dl.pdf",sep=""),height=6,width=4)
ggplot(imps,aes(var,R2))+geom_bar(stat="identity")+coord_flip()+theme_bw()
dev.off()

pdf(file=paste(pop,"_turnovers_simple_dl.pdf",sep=""),height=6,width=4)
plot_gf_turnovers(gf,imps$var[c(1:4)])
dev.off()


mm=mtrySelJack(Y=IBS.p,X=env.p,nreps=30,covariates=cbind(latlon.p.utm,admix.cov),top.pcs=tokeep,importance.cutoff = 0.1,prop.positive.cutoff = 0.5)
mm$goodvars
# dis2land did not make it
# [1] "NO2.B_range"  "N.P_range"    "NH4.B_median" "DSIGT_median" "DSIGT_range" 
# [6] "SRP.B_range"  "k490_median"  "SAL.B_range"  "dhw_max"  

save(mm,gf,file=paste(pop,"_mtrySel_dl.RData",sep=""))
load(paste(pop,"_mtrySel_post.RData",sep=""))

# boxplot of importance differences at different mtry
pop="blue"
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
                        nreps=30,top.pcs=tokeep,extra=1)
sb=ordinationJackknife(Y=IBS.p,X=cbind(env.p)[,mm$goodvars],
                       covariates=cbind(latlon.p.utm,admix.cov),
                       newX=cbind(XY,rasters.post),
                       nreps=30,top.pcs=tokeep,extra=0.1)


sum(sb.old$median.importance)
# 0.179
sum(sb$median.importance)
# 0.181

save(sb,sb.old,mm,gf,tokeep,file=paste(pop,"_rdaForest_extra01.RData",sep=""))
#load(paste(pop,"_rdaForest_extra01_post.RData",sep=""))

OFFS=gen_offset_oj(X=sb.old,Y=sb,sx=XY,sy=XY)
head(OFFS)
dim(OFFS)
# plotting the offset as a raster on the map (square coordinates):
pdf(file=paste(pop,"_offset025.pdf",sep=""),width=6,height=4)
# # simple way
plot(lonlat2raster(OFFS[,1:2],values=data.frame(cbind(offset=OFFS$offset))),col=(map.pal("viridis")))
cols=viridis(nrow(OFFS))[rank(OFFS$offset)]
plot_nice_map(OFFS[,1:2],mapdata=mapdata,margin=40000,cols=cols,pch=15,size=0.17)
dev.off()


pdf(file=paste(pop,"_turnovers_post.pdf",sep=""),height=3.5,width=3)
plot_turnover(sb,rasters.post[sb$goodrows,],"dhw_max")
plot_turnover(sb,rasters.post[sb$goodrows,],mm$goodvars[1])
plot_turnover(sb,rasters.post[sb$goodrows,],mm$goodvars[2])
plot_turnover(sb,rasters.post[sb$goodrows,],mm$goodvars[3])
plot_turnover(sb,rasters.post[sb$goodrows,],mm$goodvars[4])
plot_turnover(sb,rasters.post[sb$goodrows,],mm$goodvars[5])
dev.off()

#------ PLOTTING
# sb_old=sb
# sb=sb_new
# sb=sb_old
pdf(file=paste(pop,"_importance.pdf",sep=""),height=2.5,width=3)
important=names(sb$median.importance)
ggplot(sb$all.importances[sb$all.importances$variable %in% important,],aes(variable,importance))+
  geom_boxplot(outlier.shape = NA)+
  coord_flip()+
#  ylim(0,0.2)+
  theme_bw()
# what is the proportion of variation explained, total?
dev.off()
sum(sb$median.importance)
# Y: 0.175

pop="blue"
load(paste(pop,"_rdaForest_extra01.RData",sep=""))
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
if(is.na(meds$sum.R2[1])) {meds$sum.R2[1]=0 }

pdf(file=paste(pop,"_importance_summed.pdf",sep=""),height=2,width=2.5)
ggplot(meds,aes(predictor,sum.R2))+
  geom_bar(stat="identity")+
  #  ylim(0,0.22)+
  coord_flip()+
  #  ylim(0,0.2)+
  theme_bw()
dev.off()
sum(meds$sum.R2)
# 0.18

#-----------------
pop="blue"
load(paste(pop,"_rdaForest_extra01.RData",sep=""))
turnovers=sb$predictions.turnover
raster.vars=colnames(turnovers)

table(sb$goodrows)
xy2=XY[sb$goodrows,]
dirs=sb$predictions.direct

ras2=rasters.post[which(sb$goodrows),]
bests=names(sb$median.importance)[c(1:4)]

pdf(file=paste(pop,"_bestVars_post.pdf",sep=""),height=6,width=12)
plot(lonlat2raster(xy2,ras2[,bests]),bty="n")
dev.off()

# pa0=plot_adaptation(dirs,ras2[,bests],xy2,color.scheme="111",scal=10,cex=0.8,jitscale=0.05,rangeExp=2.5,main="unclustered RF predictions")
# 
# plot_adaptation(dirs,ras2[,bests],xy2,nclust=12,color.scheme="111",scal=10,cex=0.8,jitscale=0.05,rangeExp=2.5,merge.level=0.65, main="nclust 12 direct")
# map(database="worldHires", add=TRUE,fill=F,col="grey60")
# points(lat~lon,latlon.p,col="grey80",cex=0.8)
# 

#names(xy2)=c("lon","lat")
pdf(file=paste(pop,"_maps_extra01.pdf",sep=""),width=8,height=6)
plot_adaptation(dirs,ras2[,bests],xy2,color.scheme="000",scal=8,cex=0.8,jitscale=0.05,rangeExp=1.5)
pa=plot_adaptation(dirs,ras2[,bests],xy2,nclust=12,
                   cluster.guide = turnovers,cluster.merge.level = 0.6,
                   cluster.labels=T,color.scheme="101",
                   scal=8,cex=0.4,rangeExp=2.5,
                   jitscale=0.025,lighten=0.3)
plot_nice_map(xy2,size=0.31,mapdata=mapdata,map.on.top=T,cols=pa$colors,pch=15)
plot_nice_map(xy2,size=0.31,mapdata=mapdata,map.on.top=T,cols=pa$colors,pch=15,overlay.points = latlon[choices,2:1],pch.points = 3,cols.points ="grey40")
dev.off()

