## Phytophthora pathogen Example ----


# Data comes from a study in agriculture. The primary aim is
# the study of the spatial structure of the illness induced by
# the presence of  "Phytophthora capsici" 
# on bell peppers plants. As a second question we want to assess if
# soil umidity levels (moisture) helps in predicting the illness presence.

# data are on a regular  20x20 grid 



# first we characterize the spatial structure and then we'll add the soil umidity




# 1) How can we explore the spatial structure of these data? ----

# we need the  spdep library

require(spdep)


# read data

phytoph<-read.table("phytoph.txt",header=FALSE)	# 
phytoph.dat<-phytoph
class(phytoph)		# type of object
dim(phytoph)		#object size
colnames(phytoph)	# variables names
summary(phytoph)		# some variables information


colnames(phytoph)<-c("x","y","disease","moisture")	# better names
colnames(phytoph)

# Visualization ----
attach(phytoph)
#pdf("phytoph1.pdf")

par(mfrow=c(1,1), pty="m")	# open the graphic window
plot(x,y,type="n", lab=c(20,20,7),	xlab="X", ylab="Y")		# first define axis


# Building a data posting on pdf -----

text(x,y, labels=as.character(round(moisture,0))) # insert moisture values rounded to the nearest integer
								  
								   
for (i in 1:400){						   # here we add boxes when the illness is present 
if (disease[i]==1) 					   # 
    polygon (x[i]+c(-0.4, -0.4, 0.4, 0.4), 	   
    y[i]+c(-0.45, 0.4, 0.4, -0.45), density=0)	   
} 

title("Data posting of moisture values", cex.main=1.5)
mtext(side=3, line=0.3, "boxes highlight the phytophtora capsici presence",
	 cex=1.0)
#dev.off()


plot(x,y,type="n", lab=c(20,20,7),	xlab="X", ylab="Y")		# first define axis


# Building a data posting with moisture-----

text(x,y, labels=as.character(round(moisture,0))) # insert moisture values rounded to the nearest integer


for (i in 1:400){						   # here we color in red when the illness is present 
  if (disease[i]==1) 					   # 
    text(x[i],y[i], labels=as.character(round(moisture[i],0)),col=2)  
} 

title("Data posting of moisture values", cex.main=1.5)
mtext(side=3, line=0.3, "red numbers highlight the phytophtora capsici presence",
      cex=1.0)
# We explore the spatial structure by defining several neighborhood structure ----
# We start with a first order structure ----


# Define a spatial object

coordinates(phytoph)<-~x+y	# 
gridded(phytoph)<-T		# clarify it is a regular grid as it speeds up the computations


# Usueful functions

# knearneigh (K nearest neighbours for spatial weights)
# to build a neighborhood structure using the nearest neighbor approach
#To implement a first order neighborhood structure we set the number of neighbors to 4
#we are going to see that some problems arise at the border of the grid



# Rook (k=4) ----

phytophkn<-knearneigh(coordinates(phytoph),k=4)
names(phytophkn)					# it has 5 elements 
class(phytophkn)
## we need a list of class nb not a knn object :
# knn2nb (Neighbours list from knn object)



plot(phytophkn$x)
plot(knn2nb(phytophkn),phytophkn$x,add=T)
title("Rook neighborhood")


# Queen (k=8) ----

phytophkn<-knearneigh(coordinates(phytoph),k=8)
plot(phytophkn$x)
plot(knn2nb(phytophkn),phytophkn$x,add=T)
title("Queen neighborhood")

# To avoid the edge effect we may use the function dnearneigh that builds the graph strucuture 
#using the distance between points
# First verify the distance values among points according to the choosen structure
#

# dnearneigh (Neighbourhood contiguity by distance) ----


phytoph.dn<-dnearneigh(coordinates(phytoph),0,1)
plot(phytoph.dn,coords=coordinates(phytoph))
title("Rook neighborhood")


phytoph.dn<-dnearneigh(coordinates(phytoph),1,2)
plot(phytoph.dn,coords=coordinates(phytoph))
title("Queen neighborhood")



# 2) How can we estimate the autologistic model? ----


# Once we have the neighborhood structure we can compute the sum of neighbors
# We have to adjust for edge effects
# and eventually we estimate the pseudolikelihood maximum


# First we estimate the simple model including only a costant and the sum of neighbors ----
#We know that given the presence of spatial correlation in the data it would be incorrect to estimate a simple logistic model
.# Indeed the presence of dependence in the data is such that if we estimate the model using the logistic likelihood
#estimates are not consistent


# Geman & Graffigne ("Markov random field image models and their applications
# to computer vision", 1987) prooved that the PSL estimates are consistent for autologistic models




phymat<-nb2mat(phytoph.dn, style="B")	# we the adjacency matrix this function transforms a neighborhood list into a
							# weight matrix, using the option B we have weights =1/0
							

n<-nrow(phytoph)		# sizeof the sample: 20x20
n1i<-rep(0,n)					# number of sick neighbors
numnbr<-rep(0,n)					# number of neighbors
							# 

for (i in 1:n){					# Compute the number of neighbors 
  n1i[i]<-sum(phymat[,i]*phytoph$disease) # in each grid cell
  numnbr[i]<-sum(phymat[,i])			
}							

ni<-4*n1i/numnbr					# edge effect correction
							# 

# The edge effect happens because sites on the grid borders have a smaller number of neighbors-----
# with respect to those in center of the grid
# There are several techniques "adjusting for edge effects" the most commonly used
# is to reweight neighbors.

# For example if a site has  3 neighbors instead of the required 4 the edge correction is 4/number of neighbors
# The null model:

phy1.null<-glm(disease~ni, data=phytoph,  family=binomial)

summary(phy1.null)


# Model with moisture and first order neighborhood structure:
# AL(1) 

phy1.mois<-glm(disease~ni+moisture, data=phytoph, family=binomial)

summary(phy1.mois)

AIC(phy1.null,phy1.mois)
phy<-glm(disease~moisture, data=phytoph, family=binomial)
AIC(phy1.null,phy1.mois,phy)

#Compare with a simple logistic model
### how to compare----

##one choice is to assess the predictive ability of the model, how the model is able to detect observed presences:

## 1. get the fitted probabilities
## 2. choose a threshold for example p>0.5 imples y=1

pp=predict(phy1.null,type="response")
fit1=ifelse(pp>=0.5,1,0)

##3. compare observed and fitted
table(phytoph$disease,fit1)
##33 false positive
## 32 false negative
## 213+122 out of 400 observations correct (84.75%)

pp=predict(phy1.mois,type="response")
fit1.mois=ifelse(pp>=0.5,1,0)

##3. compare observed and fitted
table(phytoph$disease,fit1.mois)

## 27 false positive
## 45 false negative
##219+109 correct answers (82%)

#a small reduction in the prediction ability

# Now we change the neighborhood structure 
# a) second order neighborhood structure
# AL(2)

phytoph2.dn<-dnearneigh(coordinates(phytoph),0,1.5)
plot(phytoph2.dn,coords=coordinates(phytoph))
title("Second order rook neighborhood")


phymat2<-nb2mat(phytoph2.dn, style="B")	# Build the weight matrix

n2i<-rep(0,n)					# number of sick neighbors
numnbr<-rep(0,n)					# number of neighbors
							# 
for (i in 1:n){					# compute the number of sick neighbors
  n2i[i]<-sum(phymat2[,i]*phytoph$disease) # for each cell 
  numnbr[i]<-sum(phymat2[,i])			
}							

ni2<-8*n1i/numnbr					# edge effect correction
							# 

#Null-model

phy2.null<-glm(disease~ni2, data=phytoph,  family=binomial)

summary(phy2.null)


#Add moisture

phy2.mois<-glm(disease~ni2+moisture, data=phytoph, family=binomial)

summary(phy2.mois)

AIC(phy2.null,phy2.mois)
# compare:
pp=predict(phy2.null,type="response")
fit2=ifelse(pp>0.5,1,0)
table(disease,fit2)
pp=predict(phy2.mois,type="response")
fit2.mois=ifelse(pp>0.5,1,0)

table(phytoph$disease,fit2) #nulla model
#no changes wrt the first order neighborhood

table(phytoph$disease,fit2.mois)

# results visualization

par(mfrow=c(1,1), pty="m")	# 

plot(x,y,type="n", lab=c(20,20,7),	xlab="X", ylab="Y")		#set axis


# data posting 

text(x,y, labels=as.character(round(moisture,0))) # rounded moisture values
								  # 
								   
for (i in 1:400){						   # boxes where the illness is present
if (disease[i]==1) 					   #  
    polygon (x[i]+c(-0.4, -0.4, 0.4, 0.4), 	   # 
    y[i]+c(-0.45, 0.4, 0.4, -0.45), density=0)	   #
} 

title("Moisture data posting", cex.main=1.5)
mtext(side=3, line=0.3, "boxes highlight the presence of the illness",cex=1.0)




plot(x,y,type="n", lab=c(20,20,7), xlab="X", ylab="Y")
text(x,y, format(round(moisture,0)))
## 
for (i in 1:400){						   
if (fit1[i]==1) 					   
    polygon (x[i]+c(-0.4, -0.4, 0.4, 0.4), y[i]+c(-0.45, 0.4, 0.4, -0.45), density=0)	   
} 

title("Moisture data posting", cex.main=1.5)
mtext(side=3, line=0.3, "boxes highlight estimated presence",
	 cex=1.0)

text(x,y, format(round(moisture,0)),col=ifelse(disease==1,2,1))

# Compare data postings instead of tables----


# What happess if we do not consider the spatial dependence?
# 

phy.nspat<-glm(disease ~ moisture, data=phytoph, family=binomial)

pred.nspat<-fitted(phy.nspat)
pred.nspat<-ifelse (pred.nspat>=0.5,1,0)

plot(x,y,type="n", lab=c(20,20,7), xlab="X", ylab="Y")
text(x,y, format(round(moisture,0)))

for (i in 1:400){						   
if (pred.nspat[i]==1) 					   
    polygon (x[i]+c(-0.4, -0.4, 0.4, 0.4), 	   
    y[i]+c(-0.45, 0.4, 0.4, -0.45), density=0)	   
} 




title("Mositure data posting", cex.main=1.5)
mtext(side=3, line=0.3, "boxes highlight estimated presence",
	 cex=1.0)


sum(phy.nspat$y)

pre.oss<-sum(phy.nspat$y)
abs.oss<-sum(phy.nspat$y==0)

pres.pred<-sum(fitted(phy.nspat)>=0.5)
ass.pred<-sum(fitted(phy.nspat)<0.5)

fit.nspat<-ifelse(fitted(phy.nspat)>=0.5,1,0)


#phytoph<-read.table("phytoph.txt",header=FALSE)



table(phytoph$disease,fit.nspat)

#31 false positive
#84 false negative
#285 correct answers (71.2%) 


# A different package and a different model ----

require(ngspatial)

### this package estimates several models among which we have a 
## centered autologistic model
## The model that the function autologistic estimates differs from the one estimated above as it centers the 
## term of the sum of neighbors with respect to the independence expectation:
##logit(y_i)=x_i'beta+psi\sum_j (y_j -mu_j)
### where
##mu_j=1/(1+exp(-x_i'beta))
## nodes= how many cores to be used in parallel computation
## MCSE= Monte Carlo standandard error
A=adjacency.matrix(20)
  mod1<-autologistic(phytoph$disease~coordinates(phytoph)-1,phytoph,A,method="PL", control=list(nodes=2,confint="bootstrap",bootit=200))
summary(mod1)

mod2<-autologistic(phytoph$disease~coordinates(phytoph)-1,phytoph,A,method="PL", control=list(nodes=4,confint="bootstrap",bootit=200))
summary(mod2)
