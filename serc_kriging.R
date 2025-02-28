rm(list=ls())
setwd("../FL_Keys_rasters_dec2024")
  li=load("serc_2024_summaries.RData")
load("variogram_model_satPC1.RData")
li
load("keys.polygon.RData")
source("smooth.krig.R")

apply(keys.polygon,2,range)
# lon      lat
# [1,] -83.15715 24.30274
# [2,] -80.01936 25.68266

# s.pre : SERC surface variables, 1995-2008
names(b.pre)
# [1] "DEPTH_median"   "NOX.S_median"   "NO3_S_median"   "NO2.S_median"  
# [5] "NH4.S_median"   "TN.S_median"    "DIN.S_median"   "TON.S_median"  
# [9] "TP.S_median"    "SRP.S_median"   "CHLA.S_median"  "TOC.S_median"  
# [13] "SiO2.S_median"  "TURB.S_median"  "SAL.S_median"   "TEMP.S_median" 
# [17] "DO.S_median"    "TN.TP_median"   "N.P_median"     "DIN.TP_median" 
# [21] "Si.DIN_median"  "X.SAT.S_median" "DSIGT_median"   "DEPTH_range"   
# [25] "NOX.S_range"    "NO3_S_range"    "NO2.S_range"    "NH4.S_range"   
# [29] "TN.S_range"     "DIN.S_range"    "TON.S_range"    "TP.S_range"    
# [33] "SRP.S_range"    "CHLA.S_range"   "TOC.S_range"    "SiO2.S_range"  
# [37] "TURB.S_range"   "SAL.S_range"    "TEMP.S_range"   "DO.S_range"    
# [41] "TN.TP_range"    "N.P_range"      "DIN.TP_range"   "Si.DIN_range"  
# [45] "X.SAT.S_range"  "DSIGT_range" 

summary(b.post)

skrig=function(latlon,dset,polygon,vm=NULL,grid.res=1000,smooth.range=3000,verbose=FALSE){
  # smooth.range can be a vector of length 
  if(length(smooth.range)==1) { 
    sr=rep(smooth.range,ncol(dset))
  } else {
      if(length(smooth.range)!=ncol(dset)) {
        stop("Error: smooth.range should be either a single number, or a vector of length ncol(dset)")
      } else { 
        sr= smooth.range
          }
    }
  rasters=list()
  for (i in 1:ncol(dset)){
    v=names(dset)[i]
    message(v)
    sk=smooth.krig(latlon,dset,V=v,polygon=keys.polygon,variogram_model=vm,verbose=verbose,
                   grid.res=grid.res,smooth.range=sr[i])
  #  print(head(sk))
    rasters[[i]]=sk[,1]
    names(rasters)[i]=v
  }
  rasters=data.frame(do.call(cbind,rasters))
  grid=data.frame(sk[,2:3])
  return(list(rasters,grid))
}

# 
# #------- kriging s_pre
# 
 # dataset=s.pre
 # sm.r=rep(2000,ncol(dataset))
# sm.r[names(dataset)=="TURB.S_median"]=10
# sm.r[names(dataset)=="TURB.S_range"]=10
# sm.r[names(dataset)=="Si.DIN_median"]=3000
# ll=latlon.pre.s
# pdf("krig.s.5k.pdf",height=8,width=12)
# par(mfrow=c(3,4))
# skk5k.pre.s=skrig(ll,dataset,keys.polygon,smooth.range=sm.r,verbose=T)
# dev.off()
# 
# 
# #------- kriging s_post
# 
# dataset=s.post
# sm.r=rep(5000,ncol(dataset))
# sm.r[names(dataset)=="TURB.S_median"]=10
# sm.r[names(dataset)=="TURB.S_range"]=10
# sm.r[names(dataset)=="Si.DIN_median"]=3000
# ll=latlon.post.s
# pdf("krig.s.p.5k.pdf",height=8,width=12)
# par(mfrow=c(3,4))
# skk5k.post.s=skrig(ll,dataset,keys.polygon,smooth.range=sm.r,verbose=T)
# dev.off()
# 
# plot(skk5k.post.s[[1]]$SAL.S_median~skk5k.pre.s[[1]]$SAL.S_median )
# 
# # ----- comparing individual ggplots for pre and posts
# para="TURB.S_median"
# lims=range(c(skk5k.pre.s[[1]][,para],skk5k.post.s[[1]][,para]))
# XY=skk5k.pre.s[[2]]
# dd=skk5k.pre.s[[1]][,para]
# d2plot=data.frame(cbind(XY,parameter=dd))
# gg=ggplot(as.data.frame(d2plot), aes(x = lon, y = lat, color = parameter)) +
#   geom_point(size = 0.1) +
#   scale_color_viridis(limits=lims) +
#   coord_equal()+
#   ggtitle(para)+
#   theme_minimal()
# plot(gg)

#----------kriging bottom pre

dataset=b.pre
sm.r=rep(5000,ncol(dataset))
ll=latlon.pre.b
# sm.r[names(dataset)=="TURB.S_median"]=10
# sm.r[names(dataset)=="TURB.S_range"]=10
 sm.r[names(dataset)=="Si.DIN_median"]=3000
# sm.r[names(dataset)=="X.SAT_B_range"]=2000
# sm.r[names(dataset)=="NO3_B_range"]=3000
# sm.r[names(dataset)=="NO2_B_range"]=10
# sm.r[names(dataset)=="SiO2_B_range"]=2000

# sk=smooth.krig(ll,dataset,V="Si.DIN_median",polygon=keys.polygon,verbose=TRUE,
#                grid.res=1000,smooth.range=1)
# skvm=smooth.krig(ll,dataset,V="Si.DIN_median",variogram_model=variogram_model, polygon=keys.polygon,verbose=TRUE,
#                grid.res=1000,smooth.range=5000)


pdf("krig.b.pre_5k.pdf",height=8,width=12)
par(mfrow=c(3,4))
skk5k.pre.b=skrig(ll,dataset,vm=variogram_model,keys.polygon,smooth.range=sm.r,verbose=T)
dev.off()

rasters.pre=skk5k.pre.b[[1]]
latlon.pre=skk5k.pre.b[[2]]


# par(mfrow=c(1,1))
# sk=smooth.krig(ll,dataset,V="NOX.B_range",polygon=keys.polygon,verbose=TRUE,
#                grid.res=1000,smooth.range=10)
# 

#----------kriging bottom post

dataset=b.post
sm.r=rep(5000,ncol(dataset))
# sm.r[names(dataset)=="TURB.S_median"]=10
# sm.r[names(dataset)=="TURB.S_range"]=10
sm.r[names(dataset)=="Si.DIN_median"]=3000
# sm.r[names(dataset)=="X.SAT_B_range"]=2000
# sm.r[names(dataset)=="NO3_B_range"]=5500
# sm.r[names(dataset)=="NO2_B_range"]=10
# sm.r[names(dataset)=="SiO2_B_range"]=2000
ll=data.frame(latlon.post.b)
pdf("krig.b.post.5k_xy.pdf",height=8,width=12)
par(mfrow=c(3,4))
skk5kxy.post.b=skrig(ll,dataset,keys.polygon,smooth.range=sm.r,verbose=T)
dev.off()

# sk=smooth.krig(ll,dataset,V="X.SAT_B_range",polygon=keys.polygon,verbose=TRUE,
#                 grid.res=1000,smooth.range=2000)
# hist(b.post$X.SAT_B_range)


rasters.post=skk5kxy.post.b[[1]]
#rasters.pre=skk5k.pre.b[[1]]
#latlon.pre=skk5k.pre.b[[2]]
latlon.post=skk5kxy.post.b[[2]]
#identical(latlon.pre,latlon.post)
# TRUE
XY=latlon.post

save(rasters.pre,rasters.post,XY,file="SERC_krigs_bottom_prepost_5k.RData")
#save(rasters.post,XY,file="SERC_krigs_bottom_post_5k.RData")

