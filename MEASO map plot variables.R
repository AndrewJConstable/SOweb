# Script to plot variables - current and future states - for each MEASO area

library(ggplot2)
library(ggforce)
# Surface Temperature - November-December


# 


Xlim<-c(0,100)
Ylim<-c(0,25)
# all inputs are proportions of max X and Y

# depth ####

oD <- list(rect = list(xmax=0.4,ymin=0,h=1.0)
           ,XpropPerLog10 = 0.03
           ,XstartLog10 = 8)

Depth<-list(area  = 154333000
           ,depth = c(1000,2000,3000,4000)
           ,prop  = c(0.1,0.05,0.05,0.2))

# light ####

Light<-list(MidWinterPropNight = 0.7
            ,early = 20
            ,late = 60
             )

oL <- list(sun = list(x0 = 0.5, y0 = 0.8, r0 = 0.05)  # r0=radius as proportion of Xlim[2]
           ,rect = list(w = 0.05,h = 0.6)
           ,Lmin = 0
           ,Lmax = 100
           ,Ints = 10 # number of intervals
)


# temperature ####
Tearly<- 1.2
Tlate <- 5.4

oT <- list(bulb = list(x0 = 0.65, y0 = 0.2, r0 = 0.05)  # r0=radius as proportion of Xlim[2]
           ,rect = list(w = 0.05,h = 0.6)
           ,Tmin = -2
           ,Tmax = 10
           ,Ints = 11 # number of intervals
)

# mixed layer depth ####
MLD<-list(early = 100
          ,late =150)

oM <- list(wave = list(x0 = 0.8, y0 = 0.8, w = 0.1, h = 0.33, t=2) # thickness is in sine units
           ,rect = list(w = 0.05,h = 0.6)
           ,MLmin = 0
           ,MLmax = 500
           ,Ints = 10 # number of intervals
)

# sea ice ####
MLD<-list(early = 100
          ,late =150)

oM <- list(wave = list(x0 = 0.8, y0 = 0.8, w = 0.1, h = 0.33, t=2) # thickness is in sine units
           ,rect = list(w = 0.05,h = 0.6)
           ,MLmin = 0
           ,MLmax = 500
           ,Ints = 10 # number of intervals
)

########## Calculations ####

# depth polygon calculations ####

oD$box<-list(xmin = (oD$rect$xmax-(log(Depth$area)-oD$XstartLog10)*oD$XpropPerLog10)*Xlim[2]
             ,xmax = oD$rect$xmax*Xlim[2]
             ,ymin = oD$rect$ymin*Ylim[2]
             ,ymax = (oD$rect$ymin+oD$rect$h)*Ylim[2])

oD$box<-c(oD$box,list(
  ticksX = (c(seq(oD$XstartLog10,floor(log(Depth$area)),1),log(Depth$area))-oD$XstartLog10)*oD$XpropPerLog10*Xlim[2]+oD$box$xmin
  ,ticksY = (oD$rect$ymin+oD$rect$h)*Ylim[2]-seq(0,Depth$depth[length(Depth$depth)],1000)/Depth$depth[length(Depth$depth)]*oD$rect$h *Ylim[2]
  ,poly = data.frame(x=  (c(oD$XstartLog10,log(cumsum(Depth$prop)*Depth$area),oD$XstartLog10)-oD$XstartLog10)*oD$XpropPerLog10*Xlim[2]+oD$box$xmin
                     ,y= ((oD$rect$ymin+oD$rect$h)*Ylim[2]-c(0,Depth$depth,Depth$depth[length(Depth$depth)])/Depth$depth[length(Depth$depth)]*oD$rect$h*Ylim[2]))
))


# light calculations ####
oL$sun$x<-oL$sun$x0*Xlim[2]
oL$sun$y<-oL$sun$y0*Ylim[2]
oL$sun$r<-oL$sun$r0*Xlim[2]
oL$rect$xmin <- oL$sun$x-oL$rect$w*Xlim[2]/2
oL$rect$xmax <- oL$sun$x+oL$rect$w*Xlim[2]/2
oL$rect$ymax <- oL$sun$y - (oL$sun$r^2-(oL$rect$w*Xlim[2]/2)^2)^0.5
oL$rect$ymin <- oL$rect$ymax - oL$rect$h*Ylim[2]

oL$ticksY <- oL$rect$ymax-(seq(oL$Lmin,oL$Lmax,(oL$Lmax-oL$Lmin)/oL$Ints)-oL$Lmin)/(oL$Lmax-oL$Lmin)*oL$rect$h*Ylim[2]
oL$ticksX <- oL$rect$xmin+ oL$rect$w*Xlim[2]*c(0.4,0.6)

oL$poly <-data.frame(x = c(oL$rect$xmin,oL$rect$xmin,oL$rect$xmax,oL$rect$xmax)
                          ,y = c(oL$rect$ymax
                                 ,oL$rect$ymax-(Light$early-oL$Lmin)/(oL$Lmax-oL$Lmin)*oL$rect$h*Ylim[2]
                                 ,oL$rect$ymax-(Light$late-oL$Lmin)/(oL$Lmax-oL$Lmin)*oL$rect$h*Ylim[2]
                                 ,oL$rect$ymax))



# temperature thermometer calculations ####
oT$bulb$x<-oT$bulb$x0*Xlim[2]
oT$bulb$y<-oT$bulb$y0*Ylim[2]
oT$bulb$r<-oT$bulb$r0*Xlim[2]
oT$rect$xmin <- oT$bulb$x-oT$rect$w*Xlim[2]/2
oT$rect$xmax <- oT$bulb$x+oT$rect$w*Xlim[2]/2
oT$rect$ymin <- oT$bulb$y + (oT$bulb$r^2-(oT$rect$w*Xlim[2]/2)^2)^0.5
oT$rect$ymax <- oT$rect$ymin + oT$rect$h*Ylim[2]

oT$ticksY <- (seq(oT$Tmin,oT$Tmax,(oT$Tmax-oT$Tmin)/oT$Ints)-oT$Tmin)/(oT$Tmax-oT$Tmin)*oT$rect$h*Ylim[2]+oT$rect$ymin
oT$ticksX <- oT$rect$xmin+ oT$rect$w*Xlim[2]*c(0.4,0.6)

oT$tempShape <-data.frame(x = c(oT$rect$xmin,oT$rect$xmin,oT$rect$xmax,oT$rect$xmax)
                          ,y = c(oT$rect$ymin
                                 ,(Tearly-oT$Tmin)/(oT$Tmax-oT$Tmin)*oT$rect$h*Ylim[2]+oT$rect$ymin
                                 ,(Tlate-oT$Tmin)/(oT$Tmax-oT$Tmin)*oT$rect$h*Ylim[2]+oT$rect$ymin
                                 ,oT$rect$ymin))

# Mixed Layer Depth #####

oM$wave$df <- 
  {waveXt<-c(seq(0,2*pi,0.001),2*pi)
  waveXb<-c(seq(2*pi,0,-0.001),0)
  waveYt<-sin(waveXt)
  waveYb<-sin(waveXb)+oM$wave$t
  waveX<-((c(waveXt,waveXb)/(2*pi)-0.5)*oM$wave$w+oM$wave$x0)      *Xlim[2]
  waveY<-(((c(waveYt,waveYb)+1)/(2+oM$wave$t)-0.5)*oM$wave$h+oM$wave$y0)*Ylim[2]
  wave<-data.frame(x=waveX, y=waveY)
  wave}

oM<-c(oM,{
  x0<-(oM$wave$w-oM$rect$w)/2 /oM$wave$w  # proportion along wave
  xE<-(x0+oM$rect$w/oM$wave$w)
  polyX<-c(seq(x0*2*pi,xE*2*pi,0.001),xE*2*pi)
  polyY<-sin(polyX)
  polyX<-(polyX/(2*pi)-0.5)*oM$wave$w+oM$wave$x0
  polyY<-((polyY+1)/(2+oM$wave$t)-0.5)*oM$wave$h+oM$wave$y0
  rectYtop<-polyY[length(polyY)]
  rectYbottom<-polyY[length(polyY)]-oM$rect$h
  pX1<-polyX[1]
  pXe<-polyX[length(polyX)]
  pY1<-rectYtop-(MLD$early-oM$MLmin)/(oM$MLmax-oM$MLmin)*oM$rect$h
  pYe<-rectYtop-(MLD$late-oM$MLmin)/(oM$MLmax-oM$MLmin)*oM$rect$h

  boxX<-c(polyX,pXe,pX1) 
  boxY<-c(polyY,rectYbottom,rectYbottom)

  list(box = data.frame(x=boxX*Xlim[2],y=boxY*Ylim[2])
       ,poly = data.frame(x=c(pX1,polyX,pXe)*Xlim[2],y=c(pY1,polyY,pYe)*Ylim[2])
       ,ticksX = (min(polyX)+(max(polyX)-min(polyX))*c(0.4,0.6))*Xlim[2]
       ,ticksY = (rectYtop-(seq(oM$MLmin,oM$MLmax,(oM$MLmax-oM$MLmin)/oM$Ints)/(oM$MLmax-oM$MLmin)*oM$rect$h))*Ylim[2]
      )
})

# Sea ice concentration and extent
p<-ggplot() + coord_fixed(ratio=1) + xlim(Xlim) + ylim(Ylim) + xlab("")+ylab("")+
   theme(
     axis.text.x = element_blank()
     ,axis.text.y = element_blank()
     ,axis.ticks = element_blank()
     ,panel.grid.major = element_blank()
     , panel.grid.minor = element_blank()
     ,panel.background = element_rect(fill='transparent')
     ,plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
)

# depth
p<-p+geom_rect(aes(xmin=oD$box$xmin,xmax=oD$box$xmax,ymin=oD$box$ymin,ymax=oD$box$ymax),fill="white",linetype=1,colour="black",show.legend=FALSE)
p<-p+geom_polygon(data=oD$box$poly,aes(x=x,y=y),fill="brown", show.legend=FALSE)
for(i in c(1:length(oD$box$ticksX))) p<-p+geom_line(aes(x=x,y=y),data=data.frame(x=c(oD$box$ticksX[i],oD$box$ticksX[i]),y=c(oD$box$ymax-2,oD$box$ymax)),linetype=1,colour="black")
for(i in c(1:length(oD$box$ticksY))) p<-p+geom_line(aes(x=x,y=y),data=data.frame(x=c(oD$box$xmax-2,oD$box$xmax), y=c(oD$box$ticksY[i],oD$box$ticksY[i])),linetype=1,colour="black")

# light 
p<-p+ geom_circle(aes(x0 = oL$sun$x, y0 = oL$sun$y, r = oL$sun$r), fill = "yellow", show.legend=FALSE)
p<-p+geom_rect(aes(xmin=oL$rect$xmin,xmax=oL$rect$xmax,ymin=oL$rect$ymin,ymax=oL$rect$ymax),fill="white",linetype=1,colour="black",show.legend=FALSE)
p<-p+geom_polygon(data=oL$poly,aes(x=x,y=y),fill="yellow", show.legend=FALSE)
for(i in c(1:length(oL$ticksY))) p<-p+geom_line(aes(x=x,y=y),data=data.frame(x=oL$ticksX,y=c(oL$ticksY[i],oL$ticksY[i])),linetype=1,colour="black")


# thermometer 
  p<-p+ geom_circle(aes(x0 = oT$bulb$x, y0 = oT$bulb$y, r = oT$bulb$r), fill = "red", show.legend=FALSE)
  p<-p+geom_rect(aes(xmin=oT$rect$xmin,xmax=oT$rect$xmax,ymin=oT$rect$ymin,ymax=oT$rect$ymax),fill="white",linetype=1,colour="black",show.legend=FALSE)
  p<-p+geom_polygon(data=oT$tempShape,aes(x=x,y=y),fill="red", show.legend=FALSE)
    for(i in c(1:length(oT$ticksY))) p<-p+geom_line(aes(x=x,y=y),data=data.frame(x=oT$ticksX,y=c(oT$ticksY[i],oT$ticksY[i])),linetype=1,colour="black")

  # MLD
  p<-p+geom_polygon(data=oM$wave$df,aes(x=x,y=y),fill="blue",linetype=1,colour="black", show.legend=FALSE)
  p<-p+geom_polygon(data=oM$box,aes(x=x,y=y),fill="white",linetype=1,colour="black", show.legend=FALSE)
  p<-p+geom_polygon(data=oM$poly,aes(x=x,y=y),fill="blue", show.legend=FALSE)
  for(i in c(1:length(oM$ticksY))) p<-p+geom_line(aes(x=x,y=y),data=data.frame(x=oM$ticksX,y=c(oM$ticksY[i],oM$ticksY[i])),linetype=1,colour="black")
  
p

plot(oM$poly$x,oM$poly$y,type="l")

