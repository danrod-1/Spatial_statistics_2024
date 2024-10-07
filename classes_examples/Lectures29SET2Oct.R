library(geoR)
data(wolfcamp)


class(wolfcamp)
methods(class="geodata")
plot(wolfcamp)
points(wolfcamp)
summary(wolfcamp)
plot(variog(wolfcamp,option="cloud"),type="p")
plot(variog(wolfcamp,option="bin", trend= "1st", max.dist = 200),type="l")
vv<-variog(wolfcamp,option="bin", trend= "1st", max.dist = 200)
class(vv)
ww<-as.data.frame(wolfcamp)
require(tidyverse)
ggplot(ww)+
  geom_point(aes(x=X1,y=X2,size=data))+
  theme_light()
ggplot(data.frame(u=vv$u,v=vv$v))+ 
  geom_line(aes(x=u,y=v))+
  theme_light()

require(gstat)
require(sp)
?idw
data("meuse.all")
meuse<-meuse.all #save the original dataset for future  elaborations

coordinates(meuse) = ~x+y
#we need a grid where we'll interpolate observations
data("meuse.grid")
gridded(meuse.grid) = ~x+y

## rough exploration
spplot(meuse["zinc"])
meuse$lzinc<-log(meuse$zinc)
spplot(meuse["lzinc"])
spplot(meuse.grid)
#simple interpolation
x<-idw(lzinc~1,meuse, meuse.grid) #default idp=2
x1<-idw(lzinc~1,meuse, meuse.grid,idp=1.5)
x2<-idw(lzinc~1,meuse, meuse.grid,idp=3)

spplot(x["var1.pred"])
spplot(x1["var1.pred"])
spplot(x2["var1.pred"])

#surface plot can be made using library MBA and mba.surf interpolation
#mean squared error is a proper score 