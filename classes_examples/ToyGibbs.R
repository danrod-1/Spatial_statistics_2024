#####gibbs sampler -----

n=10
niter=5000
a=1
b=2
xy=data.frame(x=rep(NA,niter),y=rep(NA,niter)) 
xy$y[1] = 0.5
xy$x[1] =rbinom(1, n, xy$y[1])
xy1=data.frame(x=rep(NA,niter),y=rep(NA,niter)) 
xy1$y[1] = 0.1
xy1$x[1] =rbinom(1, n, xy1$y[1])
for(i in 2:niter){
  xy$y[i]=rbeta(1,(xy$x[(i-1)]+a), (n-xy$x[(i-1)]+ b)) 
  xy$x[i]=rbinom (1,n,xy$y[(i-1)])
  xy1$y[i]=rbeta(1,(xy1$x[(i-1)]+a), (n-xy1$x[(i-1)]+ b)) 
  xy1$x[i]=rbinom (1,n,xy1$y[(i-1)])
}
plot(xy$y,type="l")
lines(xy1$y,col=4,lty=2)
acf(xy1$y)
library(coda)
yy<-as.mcmc.list(list(as.mcmc(xy$y),as.mcmc(xy1$y)))

gelman.diag(yy)