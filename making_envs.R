# ----- forming envs for Agaricia and Porites

# ll=load("~/Dropbox/methods_in_ecogeno_2023/Agaricia_RDA_gradientForest_practice/agaricia_redo_july11_2023/Agaricia_envdata_clean_oldDepth_gebco_dec2024.RData")
# dim(env)
# a.depth=env$depth
# a.ll=latlon
# ll=load("../porites_2024/Porites_env_rasters_noSmooth_fixVar_gebco.RData")
# dim(env)
# p.depth=env$depth
# p.ll=latlon
# 
# save(p.depth,p.ll,a.depth,a.ll,file="agaricia_porites_depths_latlons.RData")

ll=load("../FL_Keys_rasters_dec2024/agaricia_porites_depths_latlons.RData")
dim(p.ll)
dim(a.ll)
# "p.depth" "p.ll"    "a.depth" "a.ll"   
ll=load("../FL_Keys_rasters_dec2024/rasters_keys_jan2025_gebcoDepth.RData")

names(rasters.pre)
# "rasters.pre"  "rasters.post" "XY"
source("~/Dropbox/keys_rdaforest_december2024/FL_Keys_rasters_dec2024/lonlat2raster.R")


sites=unique(rbind(a.ll,p.ll))

rpost=lonlat2raster(XY,rasters.post)
s.env=terra::extract(rpost,sites[,2:1])[,-1]

#------ correlations of predictors across sampled sites

s.env$DEPTH_median=NULL
s.env$DEPTH_range=NULL

sc=abs(cor(s.env))
pdf("env_cortree_feb25_completeLink.pdf",height=6,width=8)
plot(hclust(as.dist(1-sc),method="complete"),cex=0.8)
abline(h=0.1,col="red")
abline(h=0.2,col="gold")
dev.off()

s.env$TON.B_median=NULL
s.env$TON.B_range=NULL
s.env$NOX.B_median=NULL
s.env$NOX.B_range=NULL
s.env$DO.B_median=NULL
s.env$DO.B_range=NULL
s.env$Si.DIN_median=NULL
s.env$Si.DIN_range=NULL
s.env$X.SAT_B_median=NULL
s.env$DIN.B_range=NULL

sc=abs(cor(s.env))
plot(hclust(as.dist(1-sc),method="com"))
abline(h=0.1,col="red")
abline(h=0.2,col="gold")
plot(hclust(as.dist(1-sc),method="ave"))

goodenvs=names(s.env)
rpre=lonlat2raster(XY,rasters.pre[,goodenvs])
rpost=lonlat2raster(XY,rasters.post[,goodenvs])

a.env=terra::extract(rpost,a.ll[,2:1])[,-1]
p.env=terra::extract(rpost,p.ll[,2:1])[,-1]
p.env$depth=p.depth
a.env$depth=a.depth

si=unique(rbind(a.env,p.env))
sc=abs(cor(si))
plot(hclust(as.dist(1-sc),method="complete"))
abline(h=0.1,col="red")
abline(h=0.2,col="gold")


env=a.env
latlon=a.ll
save(env,latlon,file="Agaricia_env_ll_jan2025.RData")

env=p.env
latlon=p.ll
save(env,latlon,file="Porites_env_ll_jan2025.RData")



