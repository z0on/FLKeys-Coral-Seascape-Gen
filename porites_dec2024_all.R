rm(list=ls())
library(RDAforest)
library(rnaturalearth)
library(rnaturalearthdata)
library(terra)
library(viridis)

mapdata= ne_countries(scale="large",continent="north america")
#source("~/Dropbox/keys_rdaforest_december2024/FL_Keys_rasters_dec2024/smooth.krig.R")
source("~/Dropbox/keys_rdaforest_december2024/FL_Keys_rasters_dec2024/lonlat2raster.R")

ll=load("~/Dropbox/keys_rdaforest_december2024/FL_Keys_rasters_dec2024/rasters_keys_feb2025.RData")
ll=load("~/Dropbox/keys_rdaforest_december2024/FL_Keys_rasters_dec2024/Porites_env_feb2025.RData")

# admixture proportions between 3 clusters
admix=read.csv("Porites_admix.csv")[,-1]
samples=paste("por",sub(".bam","",row.names(env)),sep="")
row.names(admix)=row.names(env)=samples
row.names(latlon)=samples

# ------------ choosing subpop to work on
pop="all"
IBS.p=as.matrix(read.table(paste("Po_noclones.ibsMat",sep="")))
dimnames(IBS.p)=list(samples,samples)

#choices0=row.names(env)[which(admix[,pop]>=0.95)]
#which(samples[choices] %in% c("aa212","aa527"))
# choices=scan(paste("Po", pop,"_60",sep=""),what="c")
# length(choices)
# choices=paste("por",sub(".bam","",choices),sep="")
# # reading lineage-specific IBS, making sure its samples match to main IBS
# IBS.p=as.matrix(read.table(paste("Po", pop,".ibsMat",sep="")))
# dimnames(IBS.p)=list(choices,choices)
# latlon.p=latlon[choices,]
# admix.p=admix[choices,]
# env.p=env[choices,]
# admix.cov=admix.p[,1:3]
# admix.cov=admix.cov[,!(names(admix.cov) %in% pop)]

latlon.p=latlon
env.p=env
admix.cov=admix[,1:3]

plot(capscale(as.dist(IBS.p)~1))

mm=plot_nice_map(latlon.p,mapdata=mapdata,margin=40000,cols="cyan3",pch=3,size=2)

# converting to UTM coordinates
epsg.code=epsg.maker(mean(range(latlon.p$lon)),mean(range(latlon.p$lat)))
latlon.p.utm=latlon2UTM(latlon.p, epsg.code)
row.names(latlon.p.utm)=row.names(IBS.p)

# IBD plot
plot(as.dist(IBS.p)~dist(latlon.p.utm),pch=19,cex=0.5,col=rgb(0,0,0,0.1))

protest(capscale(dist(latlon.p.utm)~1),capscale(as.dist(IBS.p)~1))
# all: R=0.22, p=0.001
# yellow: R=0.80, p<0.001
# blue R=0.54, p<0.001
# indigo: R=0.64, p<0.001

plot(capscale(IBS.p~1+Condition(as.matrix(latlon.p.utm))),scaling=1)

# ------ hierarchical clustering tree to look for clones, wrong species collections

hc=hclust(as.dist(IBS.p),"ave")
plot(hc,cex=0.5) # clustering of samples by IBS (great to detect clones or closely related individuals)
# there are a few "low-hanging" groups of closely related individuals
# Note: to make sure that two samples are clones, their distance must be compared to the distance
# between genotyping replicates (the same sample genotyped twice)
hcut=0.2
abline(h=hcut,col="red") # this seems like a reasonable  "low-hang" threshold for calling related groups
#pheatmap::pheatmap(1-IBS.p)

#--------- retaining a single representative of each highly-related group

cuts=cutree(hc,h=hcut)
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

dim(admix.cov)

# ---- modification to ordination command when analyzing individual lineages,
#      to DISREGARD admixture with other lineages as a force shaping genetic distances
#      (conditional PCoA)

# computing conditional PCoA

ord.all=capscale(as.dist(IBS.p)~1+Condition(as.matrix(cbind(latlon.p.utm))))
plot(ord.all,scaling=1)
plot(ord.all$CA$eig/sum(ord.all$CA$eig))

# plotting simple PCA (not corrected for latlon) 
ord=capscale(as.dist(IBS.p)~1)
scc=data.frame(scores(ord,scaling=1,display="sites"))
cols=rep("gray70",nrow(admix.cov))
cols[admix.cov$yellow>=0.6]="gold"
cols[admix.cov$indigo>=0.6]="navy"
cols[admix.cov$blue>=0.6]="skyblue"

pdf("porites_all_PCA.pdf",width=4,height=4.5)
plot(scc,asp=1,col=cols,pch=16,cex=0.9,mgp=c(2.3,1,0))
dev.off()

      # scc=data.frame(scores(ord.all,scaling=1,display="sites"))
      # library(xgboost)
      # 
      # dim(rasters.post)
      # 
      # d="MDS1"
      # dtrain <- xgb.DMatrix(data = as.matrix(env.p), label = scc[,d])
      # dtest <- xgb.DMatrix(data = as.matrix(rasters.post),label=seq(1:nrow(rasters.post)))
      # 
      # # Define parameters
      # params <- list(
      #   objective = "reg:squarederror",  # Regression objective
      #   max_depth = 3,                  # Tree depth
      #   eta = 0.1,                      # Learning rate
      #   eval_metric = "rmse"            # Root Mean Squared Error
      # )
      # 
      # # Train the XGBoost model
      # xgb_model <- xgb.train(
      #   data = dtrain,
      #   params = params,
      #   nrounds = 500,                 # Number of trees
      #   watchlist = list(train = dtrain),  # Monitor performance
      # )
      # 
      # # Plot feature importance
      # importance <- xgb.importance(feature_names = colnames(as.matrix(env.p)), model = xgb_model)
      # print(importance)
      # xgb.plot.importance(importance)


# GF model 
gf=makeGF(ord.all,cbind(env.p),ntrees=1500,pcs2keep=c(1:20))
# which MDSses are predicted?
gf$result
# MDS1       MDS2       MDS8 
# 0.45428243 0.15628289 0.06150463 
# looks like we can keep the first 9 MDSes
tokeep=9

# plot importances (R2)
imps=data.frame(importance_RDAforest(gf,ord.all))
# total R2
sum(imps)
# 0.035
names(imps)="R2"
imps$var=row.names(imps)
imps$var=factor(imps$var,levels=imps$var[order(imps$R2)])
pdf(file=paste(pop,"_importances_simple_dl.pdf",sep=""),height=7,width=4)
ggplot(imps,aes(var,R2))+geom_bar(stat="identity")+coord_flip()+theme_bw()
dev.off()

pdf(file=paste(pop,"_turnovers_simple.pdf",sep=""),height=6,width=4)
plot_gf_turnovers(gf,imps$var[c(1:4)])
dev.off()

mm=mtrySelJack(Y=IBS.p,X=cbind(env.p),nreps=30,covariates=cbind(latlon.p.utm),top.pcs=9,importance.cutoff = 0.05,prop.positive.cutoff = 0.5)
mm$goodvars
# [1] "depth"        "dis2land"     "dhw_max"      "SAL.B_range"  "NH4.B_range" 
# [6] "SRP.B_median"

# "depth"        "SAL.B_range"  "dhw_max"      "NH4.B_range"  "SRP.B_median"

save(mm,gf,file=paste(pop,"_mtrySel_dl.RData",sep=""))
#load(paste(pop,"_mtrySel.RData",sep=""))

# boxplot of importance differences at different mtry
pop="all"
pdf(file=paste(pop,"_mtrySel_dl.pdf",sep=""),height=7,width=5)
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

#--------------- jackknifing (0.1) of mtry-passing variables plus latlon

#converting XY to UTM in case we need it for predictions
names(XY)=c("lon","lat")
# converting XY to UTM coordinates (distances will be less distorted)
epsg.code=epsg.maker(mean(range(XY$lon)),mean(range(XY$lat)))
XY.utm=latlon2UTM(XY, epsg.code)

sb.pre=ordinationJackknife(Y=IBS.p,X=env.p[,mm$goodvars],
                       covariates=latlon.p.utm,
                       newX=cbind(rasters.pre[,mm$goodvars],XY),
                       nreps=20,top.pcs=9,extra=0.1)
sb=ordinationJackknife(Y=IBS.p,X=env.p[,mm$goodvars],
                        covariates=latlon.p.utm,
                        newX=cbind(rasters.post,XY),
                        nreps=20,top.pcs=9,extra=0.1)
sum(sb$median.importance)
# 0.0356

save(sb,sb.pre,mm,gf,file=paste(pop,"_rdaForest_dl.RData",sep=""))

OFFS=gen_offset_oj(X=sb.pre,Y=sb,sx=XY,sy=XY)
head(OFFS)
dim(OFFS)
# plotting the offset as a raster on the map (square coordinates):
pdf(file=paste(pop,"_offset.pdf",sep=""),width=6,height=4)
# # simple way
plot(lonlat2raster(OFFS[,1:2],values=data.frame(cbind(offset=OFFS$offset))),col=(map.pal("viridis")))
cols=viridis(nrow(OFFS))[rank(OFFS$offset)]
plot_nice_map(OFFS[,1:2],mapdata=mapdata,margin=40000,cols=cols,pch=15,size=0.27)
dev.off()

#------ PLOTTING
pop="all"
load(paste(pop,"_rdaForest_dl.RData",sep=""))

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
# Y: 0.034

pdf(file=paste(pop,"_turnovers_dl.pdf",sep=""),height=3.5,width=3)
plot_turnover(sb,rasters.post[sb$goodrows,],"dhw_max")
plot_turnover(sb,rasters.post[sb$goodrows,],mm$goodvars[1])
plot_turnover(sb,rasters.post[sb$goodrows,],mm$goodvars[2])
plot_turnover(sb,rasters.post[sb$goodrows,],mm$goodvars[3])
plot_turnover(sb,rasters.post[sb$goodrows,],mm$goodvars[4])
dev.off()

pop="all"
load(paste(pop,"_rdaForest_dl.RData",sep=""))
phys=names(sb$median.importance)[grep("k490|DSIG",names(sb$median.importance))]
therm=names(sb$median.importance)[grep("TEMP|sst|dhw",names(sb$median.importance))]
chem=names(sb$median.importance)[!(names(sb$median.importance) %in% c(phys,therm,"depth"))]

Dis2land=sb$median.importance["dis2land"]
Depth=sb$median.importance["depth"]
Phys=sum(sb$median.importance[phys])
Therm=sum(sb$median.importance[therm])
Chem=sum(sb$median.importance[chem])
meds=data.frame(cbind(predictor=c("depth","dis2land","Therm","Phys","Chem")))
meds$sum.R2=c(Depth,Dis2land,Therm,Phys,Chem)
meds$predictor=factor(meds$predictor,levels=rev(c("depth","dis2land","Therm","Phys","Chem")))

pdf(file=paste(pop,"_importance_summed_dl.pdf",sep=""),height=2,width=2.5)
ggplot(meds,aes(predictor,sum.R2))+
  geom_bar(stat="identity")+
  coord_flip()+
  #  ylim(0,0.2)+
  theme_bw()
dev.off()
sum(meds$sum.R2)
#0.038
#-----------------
pop="all"
load(paste(pop,"_rdaForest.RData",sep=""))
turnovers=sb$predictions.turnover
raster.vars=colnames(turnovers)

table(sb$goodrows)
xy2=XY[sb$goodrows,]
dirs=sb$predictions.direct

ras2=rasters.post[which(sb$goodrows),]
bests=names(sb$median.importance)[c(1:4)]

pdf(file=paste(pop,"_bestVars.pdf",sep=""),height=6,width=12)
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
plot_adaptation(dirs,ras2[,bests],xy2,color.scheme="111",scal=8,cex=0.8,jitscale=0.05,rangeExp=1.5)
pa=plot_adaptation(dirs,ras2[,bests],xy2,nclust=12,
                   cluster.guide = turnovers,cluster.merge.level = 0.45, 
                   cluster.labels=T,color.scheme="111",
                   scal=15,cex=0.4,rangeExp=3,
                   jitscale=0.025,lighten=0.3)
plot_nice_map(xy2,size=0.31,mapdata=mapdata,map.on.top=T,cols=pa$colors,pch=15)
plot_nice_map(xy2,size=0.31,mapdata=mapdata,map.on.top=T,cols=pa$colors,pch=15,overlay.points = latlon[,2:1],pch.points = 3,cols.points ="grey30")
dev.off()

