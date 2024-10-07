#setwd("/Users/jona_air/Dropbox/appoggio_SA/EnglishCourse/ExtremeValues")
library(ismev)
library(evd)

rm(list=ls())
par(mfrow=c(1,1))
#Simulations from a Gaussian distribution
set.seed(19)
pdf("normal1.pdf")
curve(dnorm(x),,from=-4,to=4,ylab='normal distribution density')
x<-rnorm(100)
points(x,rep(0,100),pch=20)
abline(v=c(-3.5,3.5),col=4,lty=5)
dev.off()
data(sealevel)
tt=1912:1992
pdf("doverharwick.pdf")
plot(tt,sealevel$dover, ylim=range(sealevel,na.rm=T),xlab="years",ylab="Annual Maxima",type="b",main="Sea Level Maxima at Dover and Harwich",pch=20)
lines(tt,sealevel[,2],col=4)
points(tt,sealevel[,2],pch=20,col=4)
legend("bottomleft",c("Dover","Harwich"),pch=20,col=c(1,4),cex=0.7)
dev.off()
data(venice)
pdf("venice.pdf")
plot(as.numeric(rownames(venice)),venice[,1],pch=20,type="b",xlab="years",ylab="sealevel", main="Venice annual maximum sea level")
dev.off()
data(exchange)
tt=as.POSIXlt(rownames(exchange))
pdf("usuk.pdf")
plot(tt,exchange[,1],type="l",xlab="time",ylab="rate",main="US-UK rate of exchange")
dev.off()
data(oxford)
tt<-1901:1980
pdf("oxford1.pdf")
plot(tt,oxford,pch=20,xlab="anno",ylab="temperature (Farenheit)",type="b")
title(main="Annual maxima of Oxford's temperatures")
dev.off()
#####Gumbel distribution: graphical method
######## q-qplot
n<-length(oxford)
i<-1:n
x.q<--log(-log((i-0.5)/n))
x.o<-sort(oxford)
pdf("oxford2.pdf")
plot(x.q,x.o,xlab="theoretical quantiles",ylab="observed quantile",pch=20)
title(main="Quantile-quantile plot")
abline(a=c(83,4))
dev.off()
######################
###### Gumbel distribution: maximum likelihood
##########################
oxford.fit<-gum.fit(oxford)
gum.diag(oxford.fit)

########################################################
# GEV: MLE
########################################
oxford.gev.fit<-gev.fit(oxford)
pdf("diaggevoxford.pdf")
gev.diag(oxford.gev.fit)
dev.off()
lik.ratio<-2*(oxford.fit$nllh-oxford.gev.fit$nllh)
print(1-pchisq(lik.ratio,1))

######### GEVis better#################

zdat<-matrix(1:length(oxford),ncol=1)
oxford.trend.fit<-gev.fit(oxford,zdat,mul=1)
mle<-oxford.trend.fit$mle
se<-oxford.trend.fit$se
nllh<-oxford.trend.fit$nllh
gev.diag(oxford.trend.fit)
lik.ratio2<-2*(oxford.gev.fit$nllh-nllh)
print(1-pchisq(lik.ratio2,1)) #pvalue no evidence in favor of the trend model

########################################################
# 					GPD
#########################################################
data(rain)
pdf("rainexample.pdf")
plot(rain,type='l',xlab='day',ylab='rain',lty=3)
abline(h=20,col=2)
title(main="Daily amount of rain in South-Easth England 1914-1962")
dev.off()
fit.rain<-gpd.fit(rain, 30)
mle<-fit.rain$mle
pdf("rainMeanExcess.pdf")
mrl.plot(rain)
dev.off()
pdf("rainFitRange.pdf")
gpd.fitrange(rain, 0, 35)
dev.off()
pdf("rainDiag.pdf")
gpd.diag(fit.rain)
dev.off()

#################################################
# Non stazionariet?
################################################
trend<-matrix(1:length(rain),ncol=1)/length(rain)
fit2.rain<-gpd.fit(rain,threshold=30,ydat=trend,sigl=1,siglink=exp)
lik.ratio<-2*(fit.rain$nllh-fit2.rain$nllh)
print(1-pchisq(lik.ratio,1))
gpd.diag(fit2.rain)

