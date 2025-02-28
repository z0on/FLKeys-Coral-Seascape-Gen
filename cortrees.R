rd=abs(cor(rasters.post))
plot(hclust(as.dist(1-rd),method="ave"))
abline(h=0.1, col="red")

e0=unique(env)
ed=abs(cor(e0))
plot(hclust(as.dist(1-ed),method="ave"))
abline(h=0.1, col="red")

rr=cbind(rasters.pre,rasters.post)

rd=abs(cor(rr))
plot(hclust(as.dist(1-rd),method="ave"))
abline(h=0.1, col="red")
