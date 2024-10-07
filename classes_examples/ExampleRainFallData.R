### load data ----
load("RainData.RData")
library(gstat)
library(fields)
library(sp)
library(ggplot2)
# Bas and Cal contain the regions profiles
#check names
names(Bas)

##### we can remove the first column
Bas<-Bas[,-1]
##################### same wit Cal

Cal<-Cal[,-1]
#### plot profiles ----
plot(Bas,ylim=range(c(Bas[,2],Cal[,2])),xlim=range(c(Bas[,1],Cal[,1])),type="l")
lines(Cal)
# check the names in the data object 

names(dati)
##### select one year using a logical variable ----
w<-dati$anno==1970
## plot pluviometers in the specific year
points(dati$XUTM[w],dati$YUTM[w],pch=20)
summary(dati)

# we have elevation the let's plot it

plot(Bas,ylim=range(c(Bas[,2],Cal[,2])),xlim=range(c(Bas[,1],Cal[,1])),type="l")
lines(Cal)
text(dati$XUTM[w],dati$YUTM[w],labels=as.character(dati$Quota[w]),cex=0.6)

# in the quota32.RData object we have a fine grid for both regions
summary(Calabriaquota32)
# again the ID column can be removed 
quotacal<-Calabriaquota32[,-1]
points(quotacal[,1:2],pch=".")
summary(quotacal$quota)

# same work with the Basilicata region

quotabas<-Basquota32[,-1]
points(quotabas[,1:2],pch=".",col=4)
#### We also armonise the columnames

colnames(quotacal)<-colnames(quotabas)
quota<-data.frame(rbind(quotabas,quotacal))

points(quotabas[,1:2],pch=".")

summary(quota$Quota)

####### To explore the elevation we divide it in 4 classes

quota.cl<-cut(quota$Quota,breaks=quantile(quota$Quota),include.lowest=T)
eticq=as.character(levels(quota.cl))
quota.cl<-cut(quota$Quota,breaks=quantile(quota$Quota),include.lowest=T,labels=c(1:4))

ggplot()+geom_point(data = Bas, aes(x=XUTM,y=YUTM),size=0.8,col="lightgray")+
  geom_point(data=Cal,aes(x=XUTM,y=YUTM),size=0.8,col="lightgray")+
  geom_point(data=quota,aes(x=XUTM,y=YUTM,fill=quota.cl,col=quota.cl))+
  geom_point(data=dati[dati$anno==1970,],aes(x=XUTM,y=YUTM),col="black",shape = 15)+
  scale_color_manual("elevation", values=c(5:2))+
  scale_fill_manual("elevation", values=c(5:2))+
  theme_light()
  

#for estimation purposes it is better to change the coordinates UTM from meters to km


dati$XUTM<-dati$XUTM/1000
dati$YUTM<-dati$YUTM/1000
# select 1982
w<-anno==1982
dati.82<-dati[w,]
# attach dataset to write less ----
attach(dati.82)
###### Is a linear trend suitable for these data? ---- 
yy<-lm(totanno~XUTM+YUTM+Quota)
summary(yy)
## looks good 
trend.dati<-formula(~XUTM+YUTM+Quota)
#####  empirical variogram #######
coordinates(dati.82)<-~XUTM+YUTM
vv<-variogram(totanno~YUTM+XUTM+Quota,data = dati.82)
plot(vv,type="b")

##### Viusualize the observed surface  ----
xy<-interp(XUTM,YUTM,totanno)
surface(xy)
lines(Cal[,1]/1000,Cal[,2]/1000,lwd=2,col="white")
lines(Bas[,1]/1000,Bas[,2]/1000,lwd=2,col = "white")
############################
# estimate some theoretical variogram -----
####################

vv1.e.ols<-gstat::fit.variogram(variogram(totanno~YUTM+XUTM+Quota,data = dati.82),vgm("Exp"), fit.method = 6)
plot(vv,vv1.e.ols)

vv1.s.ols<-gstat::fit.variogram(variogram(totanno~YUTM+XUTM+Quota,data = dati.82),vgm("Sph"), fit.method = 6)
plot(vv,vv1.s.ols)

vv1.s.reml <- gstat::fit.variogram.reml(totanno~YUTM+XUTM+Quota,data = dati.82,model = vgm(psil= 400, "Sph",range = 30,nugget=500))
plot(vv,vv1.s.reml)
# the best fit seems with the exponential

# to decide which is the best variogram we use cross-validation ----

k1<-krige.cv(formula=totanno~YUTM+XUTM+Quota,dati.82,model = vv1.e.ols, nmax=40, nfold=10)
k2<-krige.cv(formula=totanno~YUTM+XUTM+Quota,dati.82,model = vv1.s.ols, nmax=40, nfold=10)
bubble(k1,"residual")
mean(k1$residual^2)
bubble(k2,"residual")
mean(k2$residual^2)
# crossvalidation confirms the exponential as best variogram

#########################
## we estimate the kriging surface now ----
#####################
quota.sp<-quota
quota.sp$XUTM<-quota$XUTM/1000
quota.sp$YUTM<-quota$YUTM/1000

coordinates(quota.sp) <- ~XUTM+YUTM

ksurf<-krige(formula=totanno~YUTM+XUTM+Quota,dati.82,model = vv1.e.ols,newdata = quota.sp)

spplot(ksurf["var1.pred"], main = "kriging interpolation")
spplot(ksurf["var1.var"], main = "kriging standard deviation")

# plots with ggplot
ksurf.d<-data.frame(XUTM=quota$XUTM/1000,YUTM=quota$YUTM/1000,pred=ksurf@data$var1.pred,var=ksurf@data$var1.var, std=sqrt(ksurf@data$var1.var))
ggplot(ksurf.d)+geom_point(aes(x=XUTM,y=YUTM, fill = pred,col=pred), size=2.2)+
  geom_point(data = Bas, aes(x=XUTM/1000, y=YUTM/1000),col="lightgray")+
  geom_point(data = Cal, aes(x=XUTM/1000, y=YUTM/1000),col="lightgray")+
  theme_light()+
  scale_fill_gradient("Total rain 1982",low = "yellow", high = "blue")+
  scale_color_gradient("Total rain 1982",low = "yellow", high = "blue")

# confidence intervals -----
ksurf.d$low<- ksurf.d$pred-1.96*ksurf.d$std
ksurf.d$up<- ksurf.d$pred+1.96*ksurf.d$std

# any "bad prediction"?
sum(ksurf.d$low<0)
# several intervals includes 0 and negative values
#to guarantee positive estimates we can use a log transformation but we have a further question
## given that the kriging model assumes Gaussian data, should we transform totanno? ----
?krigeTg
## 
shapiro.test(dati.82$totanno)#not Gaussian
transf<-MASS::boxcox(totanno~XUTM+YUTM+Quota,data = dati.82)
lambda<-transf$x[which(transf$y==max(transf$y))]#very close to zero that is very close to a log transformation
shapiro.test((dati.82$totanno^lambda-1)/lambda)#Gaussian 

# we have to understand the spatial variation of the box-cox transformed variable ----
w<-anno==1982
dati.82<-dati[w,]
dati.82$totT<-(dati.82$totanno^lambda-1)/lambda
coordinates(dati.82)<- ~XUTM+YUTM
vvT<-variogram(totT~XUTM+YUTM+Quota,data = dati.82)
plot(vvT)
#same shape as before different values
vv1.e.ols<-gstat::fit.variogram(variogram(totT~YUTM+XUTM+Quota,data = dati.82),vgm("Exp"), fit.method = 6)
plot(vvT,vv1.e.ols)

vv1.s.ols<-gstat::fit.variogram(variogram(totT~YUTM+XUTM+Quota,data = dati.82),vgm("Sph"), fit.method = 6)
plot(vvT,vv1.s.ols)
#again we choose using crossvalidation
k1<-krige.cv(formula=totT~YUTM+XUTM+Quota,dati.82,model = vv1.e.ols, nmax=40, nfold=10)
k2<-krige.cv(formula=totT~YUTM+XUTM+Quota,dati.82,model = vv1.s.ols, nmax=40, nfold=10)
bubble(k1,"residual")
mean(k1$residual^2)
bubble(k2,"residual")
mean(k2$residual^2)
#exponential works better
ksurfT<-krige(formula=totanno~YUTM+XUTM+Quota,dati.82,model = vv1.e.ols,newdata = quota.sp,lambda=lambda)
spplot(ksurfT["var1.pred"], main = "kriging interpolation")
#95% confidence interval
low<-ksurfT$var1.pred-1.96*sqrt(ksurfT$var1.var)
sum(low<0) #no "bad values"
up<-ksurfT$var1.pred+1.96*sqrt(ksurfT$var1.var)
#Using the interval proper score verify the improvement between transformed and untransformed kriging interpolation