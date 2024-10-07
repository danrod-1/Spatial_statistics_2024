#  CAR (Conditional AutoRegressive model) ----


#check this link https://cran.r-project.org/web/packages/spdep/vignettes/sids.html
#We are going to use the SIDS dataset following the tutorial by Roger Bivand available at the link above

library(sp)
library(maptools) #library to handle maps 
library(spdep)
library(spatialreg)
llCRS <- CRS("+proj=longlat +datum=NAD27")
nc <- st_read(system.file("shapes/sids.shp", package="spData")[1], quiet=TRUE)

#The shapefile format presupposes that you have three files with extensions *.shp, *.shx, and *.dbf, where the first contains the geometry data, the second the spa- tial index, and the third the attribute data. They are required to have the same name apart from the extension, and are read here using readShapeSpatial() 
class(nc)# a n sf polygon object
#pretty map of the NC counties and their centroids.We now plot some of the variables
plot(st_geometry(nc), axes=TRUE)
text(st_coordinates(st_centroid(st_geometry(nc), of_largest_polygon=TRUE)), label=nc$FIPSNO, cex=0.5)

summary(nc)
##neighbour relationships used in Cressie and Chan (1989) to the background map as a graph. The paper is available on the class page ----
gal_file <- system.file("weights/ncCR85.gal", package = "spData")[1]
ncCR85 <- read.gal(gal_file, region.id = nc$FIPSNO)
gal_file <- system.file("weights/ncCC89.gal", package = "spData")[1]
ncCC89 <- read.gal(gal_file, region.id = nc$FIPSNO)

plot(st_geometry(nc), border="grey")
plot(ncCC89, st_centroid(st_geometry(nc), of_largest_polygon), add=TRUE, col="blue")

#each component of the list is a vector with the index numbers of the neighbours of the county in question, so that the neighbours of the county with region.id of "37001" can be retreived by matching against the indices. 
summary(ncCC89)
## there are 2 counties with no neighbors
as.character(nc$NAME)[card(ncCC89) == 0]

#We analize data from 1974-1978, the grid is irregular
# Get data ----


data(nc.sids)				
names(nc.sids)
head(nc.sids)
#Build the transformed rates variable 

sids.ft <- sqrt(1000*(nc.sids$SID74/nc.sids$BIR74))+	  
  sqrt(1000*(nc.sids$SID74+1)/nc.sids$BIR74)   
#same transformation for the non white born alive
nwbir.ft <- sqrt(1000*(nc.sids$NWBIR74/nc.sids$BIR74))+	   
  sqrt(1000*(nc.sids$NWBIR74+1)/nc.sids$BIR74)		



# visualization using a choropleth map, Colors are related to the intensity of the phenomenon


par(mfrow = c(1, 1), pty="m")

sids.ord <- c(27,41,2,85,1,22,57,28,96,100, 
              53,43,69,34,7,95,11,52,48,81,21,90,64,98, 
              91,82,4,4,4,56,42,40,88,30,33,25,24,76,8, 
              73,13,59,26,16,63,55,72,6,86,87,39,66,54,
              83,60,74,65,46,78,38,36,68,32,70,67,31,99,
              5,93,29,80,17,97,20,14,51,77,47,89,94,12,
              50,61,79,92,71,10,3,58,75,45,84,15,37,9,44,
              19,62,18,49,23,35)

sids.dcol <- round(23*(1-(sids.ft- min(sids.ft))/diff(range(sids.ft)))+1)

sids.dcol <- sids.dcol[sids.ord]            

library(maps)                               

map("county", "north carolina", fill=T,col=topo.colors(24)[sids.dcol], ylim=c(33.4,36.7))

# we use a home made function image.legend ----

source("image.legend.r")

image.legend(-82.4,34.8,zlim=range(sids.ft),
             col=rev(topo.colors(24))) 

text(-82,34.95,"F-T SIDS Rates",cex=0.9)    
mtext(side=3,line=0.5,                      
      "Freeman-Tukey Square-Root Transformed\nSIDS Rates in NC (1974-1978)",cex=1.5,col=1)

# 
# Cressie & Read  adopte a distance based neighborhood structure ----
# two counties are neighbors if their centroids are less than 30 miles a part (48.28 Km).


sids <- data.frame(long=nc.sids$lon,lat=	 nc.sids$lat,sids.ft=sids.ft,nwbir.ft=nwbir.ft,nc.sids) 					
# we get the transformed rates for non white new borns and add them to the original dataset

names(sids)		
#Cressie identifies Anson county as an outlier ----

sids <- sids[sids$CNTY.ID!=2096,]				
#convert the sids data into aspatialobject
coordinates(sids) <- ~long+lat			
coords<-coordinates(sids)				

sids.nb<-dnearneigh(coords, 0, 48.28, longlat=T) #Distance based neighborhood
summary(sids.nb)
### the same as NcCC89 neighborhood

### Spatial weights ----
# we now need the matrix of spatial weights and we have to remember that some counties have zero links

sids.W.mat <- nb2mat(sids.nb,style="B", zero.policy=T)                                   
#compute the distances between counties in the neighborhood list
sids.dist <- nbdists(sids.nb,coords,longlat=T)     
#get the vector of distances
dij <- unlist(sids.dist)                           
## elements of the weighting system
term1 <- min(dij)/dij                              
#get id of neighboring counties by row
row.id <- row(sids.W.mat)[sids.W.mat==1]      
# get id of neighboring counties by columns
col.id <- col(sids.W.mat)[sids.W.mat==1]          
#numbers of new borns in 1974
numb <- sids$BIR74                               
## second term of the weigthing system
term2 <- sqrt(numb[row.id]/numb[col.id])          
##Weights
wgts <- term1*term2                               

n <- length(sids$sids.ft)                          # sample size (n=99, Anson has been removed)
#build the weights matrix
wgtmat <- matrix(0,nrow=n,ncol=n)                  

for (i in 1:length(wgts)){
  wgtmat[col.id[i],row.id[i]] <- wgts[i]                            
}
##the matrix must be converted into a weights list to be used in some functions of spdep and related packages 
sids.W <- mat2listw(wgtmat)                        
## this weights system is not symmetrical
isSymmetric(wgtmat)                                

#further more rates depends on the counties birth rates and this induces high heterogeneity
#We need to smooth the variances of the sids rates
# we can modify the weights system that is not symmetrical
#we build a diagonal matrix with the number of new borns per county 

require(stats)
D<-diag(numb,99,99)
wdgtmat<-D%*%wgtmat

isSymmetric(wdgtmat)
sids.WD <- mat2listw(wdgtmat)

#We can now use the function spautolm from spatialreg ----
#this function estimates CAR and SAR models using a generalized weighted least square approach
## from the help: "Function taking family and weights arguments for spatial autoregression model estimation by Maximum Likelihood, using dense matrix methods, not suited to large data sets with thousands of observations. With one of the sparse matrix methods, larger numbers of observations can be handled, but the interval= argument should be set. The implementation is GLS using the single spatial coefficient value, here termed lambda, found by line search using optimize to maximise the log likelihood."

#### first model: purely spatial CAR weights for the least square are given by the total number of children born in 1974 ----

sids.nullslm <- spautolm(sids.ft ~ 1, data=sids,     family="CAR", listw=sids.WD, weights=(BIR74))     
##### second model: spatial model accounting for the presence of non white new born.  ----
### weights for the least square are given by the total number of children born in 1974
sids.raceslm <- spautolm(sids.ft ~ nwbir.ft, data=sids, family="CAR", listw=sids.WD, weights=(BIR74))                                 
### third model: a not spatial model  ----
sids.raceonly <- glm(sids.ft ~ nwbir.ft,   data=sids, weights=BIR74)                        


summary(sids.nullslm)                              
summary(sids.raceslm)                              
summary(sids.raceonly)                            

## Likelihood ratio tests ----
### to assess if spatial correlation and non white born have a significant effect we perform a likelihood ratio test using LR.sarlm 

#for the race
LR.Sarlm(sids.nullslm,sids.raceslm)   
#for the space             
LR.Sarlm(sids.raceslm,sids.raceonly)               

### some residuals analysis for the spatial model ----


par(mfrow=c(1,2),pty="s")                            

hist(sids.raceslm$fit$resid,xlab="Residuals SIDS", ylab="Frequency",cex.lab=1.6,density=12, main="Residuals histogram",cex.main=1.4)

qqnorm(scale(sids.raceslm$fit$resid),cex.lab=1.5,xlab=     "Normal ",ylab="Residuals SIDS",cex.main=1.5)        
qqline(scale(sids.raceslm$fit$resid))                     



par(mfrow=c(1,1),pty="m")  

plot(sids.raceslm$fit$fitted,sids.raceslm$fit$res,  xlab="fitted",ylab="Residuals", pch=16,cex=1.3,cex.lab=1.6,  cex.main=1.6,cex.axis=1.5)
abline(h=0)                                        

### some residuals analysis for the race-only model ----


par(mfrow=c(1,2),pty="s")                            

hist(residuals(sids.raceonly),xlab="Residuals SIDS", ylab="Frequency",cex.lab=1.6,density=12, main="Residuals histogram",cex.main=1.4)

qqnorm(scale(residuals(sids.raceonly)),cex.lab=1.5,xlab=     "Normal ",ylab="Residuals SIDS",cex.main=1.5)        
qqline(scale(residuals(sids.raceonly)))                     



par(mfrow=c(1,1),pty="m")  

plot(fitted(sids.raceonly),residuals(sids.raceonly),  xlab="fitted",ylab="Residuals", pch=16,cex=1.3,cex.lab=1.6,  cex.main=1.6,cex.axis=1.5)
abline(h=0)   


# Goodness of fit for the spatial model ----

predicted <- sids.raceslm$fit$fitted               
plot(sids$sids.ft,predicted,xlim=c(1,5.2),ylim=    
       c(1,5.2),xlab="observed",     
     ylab="fitted",cex.lab=1.5, main="SIDS rates")         
abline(0,1) 
cor(sids$sids.ft,predicted)

# Goodness of fit for the race-only model ----
predicted <- fitted(sids.raceonly)              
plot(sids$sids.ft,predicted,xlim=c(1,5.2),ylim=    
       c(1,5.2),xlab="observed",     
     ylab="fitted",cex.lab=1.5, main="SIDS rates")         
abline(0,1) 
cor(sids$sids.ft,predicted)

#the spatial model  looks a little better in terms of predictions
##visualization of results  ----                              

palette(c("lightyellow","darkgoldenrod",    
          "orange","brown"))
breaks.sids <- c(-.001,2.2,3.0,3.5,7)       
sids.ygrp <- cut(sids$sids.ft,breaks.sids)  
sids.y <- rep(0,length(sids.ygrp))          
sids.y[sids.ygrp=="(-0.001,2.2]"] <- 1      
sids.y[sids.ygrp=="(2.2,3]"] <- 2           
sids.y[sids.ygrp=="(3,3.5]"] <- 3           
sids.y[sids.ygrp=="(3.5,7]"] <- 4           
sids.y <- c(sids.y[1:84],NA,sids.y[85:99])  
sids.y2 <- sids.y[sids.ord]                 
sids.fgrp <- cut(sids.raceslm$fit$          
                   signal_trend,breaks.sids)                 
sids.f <- rep(0,length(sids.ygrp))          
sids.f[sids.fgrp=="(-0.001,2.2]"] <- 1      
sids.f[sids.fgrp=="(2.2,3]"] <- 2           
sids.f[sids.fgrp=="(3,3.5]"] <- 3           
sids.f[sids.fgrp=="(3.5,7]"] <- 4           
sids.f <- c(sids.f[1:84],NA,sids.f[85:99])  
sids.f2 <- sids.f[sids.ord]                 
sids.pgrp <- cut(sids.raceslm$fit$fit,      
                 breaks.sids)                              
sids.p <- rep(0,length(sids.ygrp))          
sids.p[sids.pgrp=="(-0.001,2.2]"] <- 1      
sids.p[sids.pgrp=="(2.2,3]"] <- 2           
sids.p[sids.pgrp=="(3,3.5]"] <- 3           
sids.p[sids.pgrp=="(3.5,7]"] <- 4           
sids.p <- c(sids.p[1:84],NA,sids.p[85:99])  
sids.p2 <- sids.p[sids.ord]                 

require(maps)                                      
par(mfrow=c(3,1))                                  
map("county","north carolina",fill=T,col=sids.y2)  
title(main="SIDS rates and transformed rates")        
legend(-82.5,35.05,legend=c("< 2.2","2.2-3.0",     
                            "3.0-3.5","> 3.5"),fill=1:4,cex=0.5)             
map("county","north carolina",fill=T,col=sids.f2)  
title(main="Estimates")                        
legend(-82.5,35.05,legend=c("< 2.2","2.2-3.0",     
                            "3.0-3.5","> 3.5"),fill=1:4,cex=0.5)             
map("county","north carolina",fill=T,col=sids.p2)  
title(main="Predicted Values")                     
legend(-82.5,35.05,legend=c("< 2.2","2.2-3.0",     
                            "3.0-3.5","> 3.5"),fill=1:4,cex=0.5)             





