library(geoR)
data(wolfcamp)
set.seed(12345)
library(R2jags)
#### Data prepartion
tobeloaded=list(Z=wolfcamp$data,
X=wolfcamp$coords[,1],
Y=wolfcamp$coords[,2],
D=as.matrix(dist(wolfcamp$coords)),
N=length(wolfcamp$data),
az=2,
bz=4000,
aw=2,
bw=4000
)
partosave=c("alfa",
"betaX",
"betaY",
"tauz",
"tauw",
"phi")
n.iter=10000
n.burn=n.iter/4
n.thin=10
model.sim=jags(tobeloaded,
#inits,
parameters.to.save=partosave,
model.file="spatial.txt",
n.chains=2,
n.iter=n.iter,
n.burnin=n.burn,
n.thin=n.thin
)
traceplot(model.sim)
round(model.sim$BUGSoutput$summary,4)
model.sim<-update(model.sim,n.iter=10000,n.burnin=1000,n.thin=10)

phi.sim=model.sim$BUGSoutput$sims.matrix[,grep("phi",colnames(model.sim$BUGSoutput$sims.matrix))]
hist(phi.sim,nclass=50,prob=T,main="decay parameter")
lines(density(phi.sim),lwd=2)

hist(3/phi.sim,nclass=50,prob=T,main="practical range")
lines(density(3/phi.sim),lwd=2)

tauw.sim=model.sim$BUGSoutput$sims.matrix[,grep("tauw",colnames(model.sim$BUGSoutput$sims.matrix))]

hist(tauw.sim,nclass=50,prob=T,main="process precision")
lines(density(tauw.sim),lwd=2)

hist(1/tauw.sim,nclass=50,prob=T,main="process varianve")
lines(density(1/tauw.sim),lwd=2)


tauz.sim=model.sim$BUGSoutput$sims.matrix[,grep("tauz",colnames(model.sim$BUGSoutput$sims.matrix))]

hist(tauz.sim,nclass=50,prob=T,main="measurement precision")
lines(density(tauz.sim),lwd=2)

hist(1/tauz.sim,nclass=50,prob=T,main="nugget")
lines(density(1/tauz.sim),lwd=2)


