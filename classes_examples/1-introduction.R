#####################################
##																###
##		Graphical exploration						###
##																###
#####################################
#library(geoR)
library(gstat)
library(spdep)
library(scatterplot3d)
library(MBA)
load("datiAll.RData")
############# scatterplot of continuous data ----

#pdf("../introductionAndGeostat/Images/scatter3d.pdf")
par(mfrow=c(1,1),pty='s',cex.lab=0.7,cex.main=1.0)
scatterplot3d(wolf$x,wolf$y,wolfcamp$piezo, type = "h", highlight.3d = TRUE,pch = 20,angle=40, 
main = "Piezometric head measurements \n at the Wolfcamp Aquifer (Texas, USA)",xlab='x',ylab='y',zlab='Piezometric head' )
#dev.off()

#################### exploratory analysis of georeferenced data ----
# We now create a function to draw 4 plots 
geoplot<-function(data){
  #data must be organized in a dataframe with first two columns = coordinates and 
  #thirs column= data
  par(mfrow=c(2,2))
  a<-cut(data[,3], breaks=quantile(data[,3]),include.lowest = T,labels = c(2,3,4,5))
  pp<-ifelse(a=="2",20,ifelse(a=="3",15,ifelse(a=="4", 18,19)))
  plot(data[,1],data[,2],xlab="x",ylab="y",col=as.numeric(as.character(a)),pch=pp)
  plot(data[,1],data[,3],xlab="x",ylab="value",pch=20)
  plot(data[,3],data[,2],ylab="y",xlab="value",pch=20)
  hist(data[,3],prob=T,main="histogram and density line")
  lines(density(data[,3]))
}
#dev.off()
geoplot(wolf)

## same plot using ggplot and gridExtra
geoplot.gg<-function(data){
  #data must be organized in a dataframe with first two columns = coordinates (x,y) and 
  #thirs column= data named value
  colnames(data)[1:3]<-c("x","y","value")
  require(ggplot2)
  require(gridExtra)
  
  val<-cut(data[,3], breaks=quantile(data[,3]),include.lowest = T,labels = c(2,3,4,5))
  #pp<-ifelse(a=="2",20,ifelse(a=="3",15,ifelse(a=="4", 18,19)))
  grid.arrange(
    ggplot(data)+geom_point(aes(x=x,y=y,fill=val,shape=val,col=val),show.legend = F)+theme_light()+xlab("x")+ylab("y"),
    ggplot(data)+geom_point(aes(x=x,y=value),show.legend = F)+theme_light()+xlab("x")+ylab("value"),
    ggplot(data)+geom_point(aes(y=x,x=value),show.legend = F)+theme_light()+ylab("y")+xlab("value"),
    ggplot(data, aes(x=x)) +
      geom_histogram(aes(y = ..density..), binwidth=density(data$x)$bw, fill=4,col=4) +
      geom_density(fill="red", alpha = 0.3)+
      theme_light()
  )
}
geoplot.gg(wolf)
#####Lattice data----
####															
library(maptools)
data(nc.sids)
######### example using SIDS data
# readShapePoly upload a raster representation of a geograohical map

sids.shp <- readShapePoly(system.file("shapes/sids.shp",package="maptools"))
sids <- sids.shp$att.data
code<-findInterval(nc.sids$SID74,trunc(seq(0,max(nc.sids$SID74),length=7)),rightmost.closed=TRUE)
sids.break <- as.ordered(cut(nc.sids$SID74,breaks=trunc(seq(0,max(nc.sids$SID74),length=7)),include.lowest=TRUE))

cols<-rainbow(7)
#pdf("../introductionAndGeostat/Images/sidsmap.pdf")
par(mfrow=c(1,1),pty='s',cex.lab=0.7,cex.main=1.0)
plot(sids.shp, col=cols[code],border="black") 
legend(c(-84,-80),c(32,34), legend=paste("counts", levels(sids.break)),  fill=cols, bty="n",cex=0.7)
 title("Sudden infant deaths \n in North Carolina for 1974-78 ")  
 #dev.off()
 
 ########################### Point processes ----
 ## not environmental data just an example to make life easier 
 library(ggmap)
 
 ######## the Houston crime data
 data(crime)
 ### subset the data and get only murder
 murder <- subset(crime, offense == "murder")
 # plot on a google map
 # pdf("../introductionAndGeostat/Images/murder.pdf")
  qmplot(lon, lat, data = murder, colour = I('red'), size = I(3), darken = .3,source="stamen")
  #dev.off()
    # pdf("../introductionAndGeostat/Images/houstonmap.pdf")
  # Houston=get_map("Houston",zoom=12)
  # ggmap(Houston)
# dev.off()

qmplot(lon, lat, data = murder, colour = I('red'), size = I(3), darken = .3,source="stamen")+geom_segment(aes(x=-95.45,xend=-95.45,y=29.7,yend=29.85))+geom_segment(aes(x=-95.40,xend=-95.40,y=29.7,yend=29.85))+geom_segment(aes(x=-95.35,xend=-95.35,y=29.7,yend=29.85))+geom_segment(aes(x=-95.30,xend=-95.30,y=29.7,yend=29.85))+geom_segment(aes(x=-95.45,xend=-95.30,y=29.7,yend=29.7))+geom_segment(aes(x=-95.45,xend=-95.30,y=29.75,yend=29.75))+geom_segment(aes(x=-95.45,xend=-95.30,y=29.8,yend=29.8))+geom_segment(aes(x=-95.45,xend=-95.30,y=29.85,yend=29.85))

#dev.off()
