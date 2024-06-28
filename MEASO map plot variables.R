# Script to plot variables - current and future states - for each MEASO area
# by Andrew J. Constable
# April 2024

# Functions to plot graphic summarising depth and for two times (early, late) the following variables:
#     sea ice, light, temperature, MLD

# 1.0 Notes ####


# 2.0 Setup for trials ####

#           2.1 Libraries ####

library(ggplot2)
library(ggforce)
library(terra)

#   2.2 Layout of one plot ####
#       2.2.1 Object (o) dimensions ####

# Depth
oD <- list(rect = list(xmax=0.5,ymin=0,h=0.75) # depth is plotted to the left of xmax, the width of x is dependent on log10 area
           ,unit         = 1E6   # if NA then log10 transform
           ,XpropPerUnit = 0.04 # proportion of X range of box per unit (on log10 scale or per unit)
           ,Xstart       = 0)    # if  then no log transform otherwise units of area on log10 scale

# Light
oL <- list(sun = list(x0 = 0.6, y0 = 0.8, r0 = 0.05, rayLen = 0.03, rayN =12)  # x0,y0 = centre of circle as proportions of ranges
                                                                               # r0=radius as proportion of Xlim[2]
                                                                               # rayLen as proportion of r0; rayN = number of rays 
           ,rect = list(w = 0.05,h = 0.6) # proportions of X and Y ranges 
           ,Lmin = 200                      # plotting values = minimum light
           ,Lmax = 280                    # maximum light
           ,Ints = 8                     # number of intervals
           ) # light

# Temperature
oT <- list(bulb = list(x0 = 0.75, y0 = 0.2, r0 = 0.05)  # as for light
           ,rect = list(w = 0.05,h = 0.6) # proportions of X and Y ranges
           ,Tmin = -2                     # minimum temperature
           ,Tmax = 14                     # maximum temperature
           ,Ints = 8                     # number of intervals
           ) # temp

# Mixed Layer Depth
oM <- list(wave = list(x0 = 0.9, y0 = 0.8, w = 0.1, h = 0.33, t=2) # as for light - thickness (t) is in sine units
           ,rect = list(w = 0.05,h = 0.6) # proportions of X and Y ranges
           ,MLmin = 30                     # minimum MLD
           ,MLmax = 60                   # maximum MLD
           ,Ints = 3                     # number of intervals
           ) # MLD

# Sea ice (on top of depth polygon)
oS <- list(rect = list(xmax = oD$rect$xmax,ymin=0.8,h=0.2)
           ,unit            = oD$unit   # if NA then log10 transform
           ,XpropPerUnit    = oD$XpropPerUnit
           ,Xstart          = oD$Xstart)

#   2.3 Trial Data ####
#       2.3.1 Plot Summary rectangle ####

Xlim<-c(0,100) # Layout dimensions are proportions of these ranges (keep minima to 0)
Ylim<-c(0,25)

#       2.3.2 Input data ####

# GEBCO data 2023

gebco23<-rast("/Users/acon/Desktop/_w/_d/GEBCO_10_Apr_2024_41650c45d0e9/GEBCO23 Southern Ocean.tif")
names(gebco23)<-"height"

# MEASO areas
rootMEASOshp<-"/Users/acon/Desktop/_w/_d/Shapefiles/MEASO polygons/"
fMEASOshp<-"MEASO_polygons.shp"

MEASOshp<-vect(paste0(rootMEASOshp,fMEASOshp))
MEASOshp<-project(MEASOshp,gebco23)

#   1     2     3     4     5     6     7     8     9     10    11    12    13    14    15
MEASOnames<-c("AOA", "AON", "AOS", "CIA", "CIN", "CIS", "EIA", "EIN", "EIS", "EPA", "EPN", "EPS", "WPA", "WPN", "WPS")
dDepth<-readRDS(file="/Users/acon/Desktop/_w/_r/SOweb/Input\ Data\ Processing/DepthProfile.rds")
MEASOsummary<-readRDS("/Users/acon/Desktop/_w/_r/SOweb/Input\ Data\ Processing/MEASOareaSummary.rds")

MEASOareaLimits<-list(minimum  = min(unlist(lapply(dDepth,function(d){d$TotalArea*min(d$prop[d$prop>0])})),na.rm=TRUE)
                      ,maximum = max(unlist(lapply(dDepth,function(d){d$TotalArea})),na.rm=TRUE))
MEASOareaLimits

MEASOlimits<-lapply(MEASOsummary,function(s){
  return(c(min(s[,ncol(s)],na.rm=TRUE),max(s[,ncol(s)],na.rm=TRUE)))
})
MEASOlimits # check limits in object data are correct

MEASOareaSummary<-lapply(c(1:15),function(m,n,d,s){
  return(list(Depth = dDepth[[m]]
             ,Light = list(early = s[["light"]][s[["light"]][,"measo"]==m,"early"]
                          ,late = s[["light"]][s[["light"]][,"measo"]==m,"late"]
                          ) # end Light
             ,Temp  = list(early = s[["temp"]][s[["temp"]][,"measo"]==m,"early"]
                          ,late = s[["temp"]][s[["temp"]][,"measo"]==m,"late"]
                          ) # end Temp
             ,MLD   = list(early = s[["MLD"]][s[["MLD"]][,"measo"]==m,"early"]
                          ,late = s[["MLD"]][s[["MLD"]][,"measo"]==m,"late"]
                          ) # end MLD
             ,SeaIce = list(MEASOarea  = dDepth[[m]]$TotalArea
                           ,conc95 = list(early = s[["SeaIce"]][s[["SeaIce"]][,"measo"]==m,"early"], late = s[["SeaIce"]][s[["SeaIce"]][,"measo"]==m,"late"])
                           ,conc75 = list(early = s[["SeaIce"]][s[["SeaIce"]][,"measo"]==m,"early"], late = s[["SeaIce"]][s[["SeaIce"]][,"measo"]==m,"late"])
                           ,conc25 = list(early = s[["SeaIce"]][s[["SeaIce"]][,"measo"]==m,"early"], late = s[["SeaIce"]][s[["SeaIce"]][,"measo"]==m,"late"])
                            ) # end SeaIce
               )) # end list and return
},MEASOnames,dDepth,MEASOsummary)


#       2.3.3 Trial plot data ####

doTicks<-list(depth=FALSE, light=FALSE,temp=FALSE, MLD=FALSE, ice=FALSE)
Colours<-list(depth = list(bg = "white", main = "black")
              ,light = list(bg = "white", main = "black")
              ,temp  = list(bg = "white", main = "black")
              ,MLD   = list(bg = "white", main = "black")
              ,ice   = list(bg = "white", c95 = "white", c75 = "grey", c25 = "black")
              )# end colours

MEASOcolRGBA<-matrix(c(2,    226, 20, 20, 1
                       ,3,    242,113,113, 1
                       ,1,    249,188,188, 1
                       ,5,    236,133, 69, 1
                       ,6,    243,177,136, 1                 
                       ,4,    249,215,194, 1
                       ,8,    135,141,199, 1                   
                       ,9,    176,180,217, 1
                       ,7,    209,211,233, 1
                       ,11,    75,125,126, 1
                       ,12,    95,158,160, 1
                       ,10,   160,197,199, 1
                       ,14,   132, 57,116, 1
                       ,15,   198,122,181, 1
                       ,13,   225,184,216, 1),ncol=5,byrow=TRUE)
dimnames(MEASOcolRGBA)[[2]]<-c("Area","R","G","B", "Alpha")
mc<-MEASOcolRGBA[order(MEASOcolRGBA[,"Area"]),]
MEASOcols<-rgb(red=mc[,"R"],green=mc[,"G"],blue=mc[,"B"],names=MEASOnames,maxColorValue = 255)

# test colours
library(colorspace)
demoplot(MEASOcols,type="bar")
detach(package:colorspace,unload=TRUE) # logical is for detaching namespace as well


# 3.0 Functions ####
plotMEASOareaSummary<-function(dIn,oD,oL,oT,oM,oS,Xlim,Ylim,doTicks,Colours,ContShelf=TRUE){ # ContShelf is to plot seabed area from left, else island/bank

  # breakup input data
  Depth  <- dIn$Depth
  Light  <- dIn$Light
  Temp   <- dIn$Temp
  MLD    <- dIn$MLD
  SeaIce <- dIn$SeaIce
  
  
# depth polygon calculations ####
Xrange<- if(is.na(oD$unit)) (log10(Depth$TotalArea)-oD$Xstart) else (Depth$TotalArea-oD$Xstart)/oD$unit
oD$box<-list(xmin = (oD$rect$xmax-Xrange*oD$XpropPerUnit)*Xlim[2]
             ,xmax = oD$rect$xmax*Xlim[2]
             ,ymin = oD$rect$ymin*Ylim[2]
             ,ymax = (oD$rect$ymin+oD$rect$h)*Ylim[2])

  # if ContShelf=TRUE then plot bottom profile from left, else split profile in centre to represent islands/banks
  # also need to scale bottom profile to max height when less than 0
  useProp<-Depth$prop>0
  X0    <- if(is.na(oD$unit)) 10^oD$Xstart else oD$Xstart
  
  dX<-c(X0,Depth$prop[useProp]*Depth$TotalArea)
  OriginDepth<-if(Depth$MaxHeight>=0) 0 else Depth$MaxHeight 
  dY<-c(OriginDepth,Depth$depth[useProp])   

  dX<-c(dX,dX[1]); dY<-c(dY,dY[length(dY)]) # returning the line to close the polygon

  if(is.na(oD$unit)){  
  dX <-(log10(dX)-oD$Xstart)*oD$XpropPerUnit*Xlim[2]+oD$box$xmin
  } else {
  dX <-  (dX-oD$Xstart)/oD$unit*oD$XpropPerUnit*Xlim[2]+oD$box$xmin
  }

  dY<-(oD$rect$ymin+oD$rect$h)*Ylim[2]-dY/Depth$depth[length(Depth$depth)]*oD$rect$h*Ylim[2]
  
xTickLast<- if(is.na(oD$unit)) floor(log10(Depth$TotalArea)-oD$Xstart) else floor((Depth$TotalArea-oD$Xstart)/oD$unit)
xTickEnd <- if(is.na(oD$unit)) (log10(Depth$TotalArea)-oD$Xstart) else (Depth$TotalArea--oD$Xstart)/oD$unit

oD$box<-c(oD$box,list(
  ticksX =  c(seq(0,xTickLast,1),xTickEnd)*oD$XpropPerUnit*Xlim[2]+oD$box$xmin # need to correct if not continental shelf
  ,ticksY = (oD$rect$ymin+oD$rect$h)*Ylim[2]-seq(0,Depth$depth[length(Depth$depth)],-1000)/(-Depth$depth[length(Depth$depth)])*oD$rect$h *Ylim[2]
  ,poly = data.frame(x= dX, y= dY)
))

# light calculations ####
oL$sun$x<-oL$sun$x0*Xlim[2]
oL$sun$y<-oL$sun$y0*Ylim[2]
oL$sun$r<-oL$sun$r0*Xlim[2]

thetaInt<-2*pi/oL$sun$rayN
theta<-seq(0,2*pi,thetaInt)

oL$sun$poly<-as.data.frame(matrix(unlist(sapply(theta,function(t,sun,L,tI){
  r<-c((sun$r-L),sun$r)
  theta<-c(t,t+tI/2)
  x<-cos(theta)*r+sun$x
  y<-sin(theta)*r+sun$y
  return(rbind(theta,x,y))
   },oL$sun,oL$sun$rayLen*Xlim[2],thetaInt)),ncol=3,byrow=TRUE))
names(oL$sun$poly)<-c("theta","x","y")

# background box to hide rays
oL$rectBG$xmin <- oL$sun$x-oL$rect$w*Xlim[2]/2
oL$rectBG$xmax <- oL$sun$x+oL$rect$w*Xlim[2]/2
oL$rectBG$ymin <- oL$rect$ymax-oL$sun$rayLen*Xlim[2] - oL$rect$h*Ylim[2]

# data box
oL$rect$xmin <- oL$sun$x-oL$rect$w*Xlim[2]/2
oL$rect$xmax <- oL$sun$x+oL$rect$w*Xlim[2]/2
oL$rect$ymax <- oL$sun$y - (oL$sun$r^2-(oL$rect$w*Xlim[2]/2)^2)^0.5
oL$rect$ymin <- oL$rect$ymax - oL$rect$h*Ylim[2]
oL$rect$ymaxBG <- oL$sun$y - ((oL$sun$r-oL$sun$rayLen*Xlim[2])^2-(oL$rect$w*Xlim[2]/2)^2)^0.5


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
                                 ,(Temp$early-oT$Tmin)/(oT$Tmax-oT$Tmin)*oT$rect$h*Ylim[2]+oT$rect$ymin
                                 ,(Temp$late-oT$Tmin)/(oT$Tmax-oT$Tmin)*oT$rect$h*Ylim[2]+oT$rect$ymin
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
       ,polybase = data.frame(x=c(polyX,pX1),y=c(polyY,rectYtop))
       ,poly = data.frame(x=c(pX1,polyX,pXe)*Xlim[2],y=c(pY1,polyY,pYe)*Ylim[2])
       ,ticksX = (min(polyX)+(max(polyX)-min(polyX))*c(0.4,0.6))*Xlim[2]
       ,ticksY = (rectYtop-(seq(oM$MLmin,oM$MLmax,(oM$MLmax-oM$MLmin)/oM$Ints)/(oM$MLmax-oM$MLmin)*oM$rect$h))*Ylim[2]
      )
})

# Sea ice concentration and extent (plotted over depth) ####
sXrange<- if(is.na(oS$unit)) (log10(SeaIce$MEASOarea)-oS$Xstart) else (SeaIce$MEASOarea-oS$Xstart)/oS$unit

oS$box<-list(xmin = (oS$rect$xmax-sXrange*oS$XpropPerUnit)*Xlim[2]
             ,xmax = oS$rect$xmax*Xlim[2]
             ,ymin = oS$rect$ymin*Ylim[2]
             ,ymax = (oS$rect$ymin+oS$rect$h)*Ylim[2])


if(is.na(oS$unit)){  
  s95<-list(early = log10(SeaIce$conc95$early/100*SeaIce$MEASOarea)-oS$Xstart
            ,late = log10(SeaIce$conc95$late/100*SeaIce$MEASOarea)-oS$Xstart)
  s75<-list(early = log10(SeaIce$conc75$early/100*SeaIce$MEASOarea)-oS$Xstart
            ,late = log10(SeaIce$conc75$late/100*SeaIce$MEASOarea)-oS$Xstart)
  s25<-list(early = log10(SeaIce$conc25$early/100*SeaIce$MEASOarea)-oS$Xstart
            ,late = log10(SeaIce$conc25$late/100*SeaIce$MEASOarea)-oS$Xstart)
} else {
  s95<-list(early = (SeaIce$conc95$early/100*SeaIce$MEASOarea-oS$Xstart)/oS$unit
            ,late = (SeaIce$conc95$late/100*SeaIce$MEASOarea-oS$Xstart)/oS$unit)
  s75<-list(early = (SeaIce$conc75$early/100*SeaIce$MEASOarea-oS$Xstart)/oS$unit
            ,late = (SeaIce$conc75$late/100*SeaIce$MEASOarea-oS$Xstart)/oS$unit)
  s25<-list(early = (SeaIce$conc25$early/100*SeaIce$MEASOarea-oS$Xstart)/oS$unit
            ,late = (SeaIce$conc25$late/100*SeaIce$MEASOarea-oS$Xstart)/oS$unit)
}

oS$box<-c(oS$box,list(
   poly95 = data.frame(x=c(oS$box$xmin,oS$box$xmin
                           ,oS$box$xmin+s95$early*oS$XpropPerUnit*Xlim[2]
                           ,oS$box$xmin+s95$late*oS$XpropPerUnit*Xlim[2]
                           )
                       ,y=c(oS$box$ymin,oS$box$ymax,oS$box$ymax,oS$box$ymin))
  ,poly75 = data.frame(x=c(oS$box$xmin,oS$box$xmin
                           ,oS$box$xmin+s75$early*oS$XpropPerUnit*Xlim[2]
                           ,oS$box$xmin+s75$late*oS$XpropPerUnit*Xlim[2]
                           )
                       ,y=c(oS$box$ymin,oS$box$ymax,oS$box$ymax,oS$box$ymin))
  ,poly25 = data.frame(x=c(oS$box$xmin,oS$box$xmin
                           ,oS$box$xmin+s25$early*oS$XpropPerUnit*Xlim[2]
                           ,oS$box$xmin+s25$late*oS$XpropPerUnit*Xlim[2]
                           )
                       ,y=c(oS$box$ymin,oS$box$ymax,oS$box$ymax,oS$box$ymin))
))

# Plot ####

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
p<-p+geom_rect(aes(xmin=oD$box$xmin,xmax=oD$box$xmax,ymin=oD$box$ymin,ymax=oD$box$ymax),fill=Colours$depth$bg,linetype=1,colour="black",show.legend=FALSE)
p<-p+geom_polygon(data=oD$box$poly,aes(x=x,y=y),fill=Colours$depth$main, show.legend=FALSE)
if(doTicks$depth) for(i in c(1:length(oD$box$ticksX))) p<-p+geom_line(aes(x=x,y=y),data=data.frame(x=c(oD$box$ticksX[i],oD$box$ticksX[i]),y=c(oD$box$ymax-2,oD$box$ymax)),linetype=1,colour="black")
if(doTicks$depth) for(i in c(1:length(oD$box$ticksY))) p<-p+geom_line(aes(x=x,y=y),data=data.frame(x=c(oD$box$xmax-2,oD$box$xmax), y=c(oD$box$ticksY[i],oD$box$ticksY[i])),linetype=1,colour="black")

# light 
#p<-p+ geom_circle(aes(x0 = oL$sun$x, y0 = oL$sun$y, r = oL$sun$r), fill = Colours$light$main, show.legend=FALSE)
p<-p+geom_polygon(data=oL$sun$poly,aes(x=x,y=y),fill=Colours$light$main,linetype=1,colour="black", show.legend=FALSE)
#p<-p+geom_rect(aes(xmin=oL$rect$xmin,xmax=oL$rect$xmax,ymin=oL$rect$ymax,ymax=oL$rect$ymaxBG),fill=Colours$light$main,show.legend=FALSE)
p<-p+geom_rect(aes(xmin=oL$rect$xmin,xmax=oL$rect$xmax,ymin=oL$rect$ymin,ymax=oL$rect$ymax),fill=Colours$light$bg,linetype=1,colour="black",show.legend=FALSE)
p<-p+geom_polygon(data=oL$poly,aes(x=x,y=y),fill=Colours$light$main, show.legend=FALSE)
if(doTicks$light) for(i in c(1:length(oL$ticksY))) p<-p+geom_line(aes(x=x,y=y),data=data.frame(x=oL$ticksX,y=c(oL$ticksY[i],oL$ticksY[i])),linetype=1,colour="black")


# thermometer 
  p<-p+ geom_circle(aes(x0 = oT$bulb$x, y0 = oT$bulb$y, r = oT$bulb$r), fill = Colours$temp$main, show.legend=FALSE)
  p<-p+geom_rect(aes(xmin=oT$rect$xmin,xmax=oT$rect$xmax,ymin=oT$rect$ymin,ymax=oT$rect$ymax),fill=Colours$temp$bg,linetype=1,colour="black",show.legend=FALSE)
  p<-p+geom_polygon(data=oT$tempShape,aes(x=x,y=y),fill=Colours$temp$main, show.legend=FALSE)
  if(doTicks$temp)   for(i in c(1:length(oT$ticksY))) p<-p+geom_line(aes(x=x,y=y),data=data.frame(x=oT$ticksX,y=c(oT$ticksY[i],oT$ticksY[i])),linetype=1,colour="black")

  # MLD
  p<-p+geom_polygon(data=oM$wave$df,aes(x=x,y=y),fill=Colours$MLD$main,linetype=1,colour="black", show.legend=FALSE)
  p<-p+geom_polygon(data=oM$box,aes(x=x,y=y),fill=Colours$MLD$bg,linetype=1,colour="black", show.legend=FALSE)
  p<-p+geom_polygon(data=oM$poly,aes(x=x,y=y),fill=Colours$MLD$main, show.legend=FALSE)
  p<-p+geom_polygon(data=oM$polybase,aes(x=x,y=y),fill=Colours$MLD$main,linetype=1,colour="black", show.legend=FALSE)
  
  if(doTicks$MLD) for(i in c(1:length(oM$ticksY))) p<-p+geom_line(aes(x=x,y=y),data=data.frame(x=oM$ticksX,y=c(oM$ticksY[i],oM$ticksY[i])),linetype=1,colour="black")

  if(sum(oS$box$poly95$x>oS$box$xmin)>0){
  # Sea ice (data are for open water)  
  # p<-p+geom_rect(aes(xmin=oS$box$xmin,xmax=oS$box$xmax,ymin=oS$box$ymin,ymax=oS$box$ymax),fill=Colours$ice$bg,linetype=1,colour="black",show.legend=FALSE)
  p<-p+geom_polygon(data=oS$box$poly95,aes(x=x,y=y),fill=Colours$ice$c95,linetype=1,colour="black", show.legend=FALSE)
  p<-p+geom_polygon(data=oS$box$poly75,aes(x=x,y=y),fill=Colours$ice$c75, show.legend=FALSE)
  p<-p+geom_polygon(data=oS$box$poly25,aes(x=x,y=y),fill=Colours$ice$c25, show.legend=FALSE)
  if(doTicks$ice) for(i in c(1:length(oS$box$ticksX))) p<-p+geom_line(aes(x=x,y=y),data=data.frame(x=c(oS$box$ticksX[i],oS$box$ticksX[i]),y=c(oS$box$ymax-2,oS$box$ymax)),linetype=1,colour="black")
  }
  

return(p)
}

fnJoinWParea<-function(gt_m180,lt_m180,r,shp){  # input polygon numbers
  s1<-r[r[,"geom"]==gt_m180,]
  s2<-r[r[,"geom"]==lt_m180,]
  #plot(vect(s1,type="polygons",crs=crs(shp)))
  s1180y<-s1[s1[,"x"]<(-179.999),"y"]
  s2180y<-s2[s2[,"x"]>(-180.001),"y"]
  
  s1Mod<-s1[which(!(s1[,"x"]<(-179.999) & !(s1[,"y"]==min(s1180y) | s1[,"y"]==max(s1180y)))),]
  s2Mod<-s2[which(!(s2[,"x"]>(-180.001) & !(s2[,"y"]==min(s2180y) | s2[,"y"]==max(s2180y)))),]
  
  s1Mod<-s1Mod[-1,]
  
  p180<-which(s2Mod[,"x"]>(-180.001))
  sFull<-rbind(s2Mod[c(1:p180[1]),],s1Mod,s2Mod[-c(1:p180[1]),])
  sFull[,"part"]<-1
  sFull[,"geom"]<-1
  
  return(vect(sFull,type="polygons",crs=crs(shp)))
}

plotMEASObaseMap<-function(shp){
  
  # for each of the west pacific areas (shapes 13:15), part 2 of the shape is to be
  # at each end of the plot.  Those on the left are to be joined to part 1 as a single
  # polygon and those on the right would be separate polygons
  
  r0 <- geom(shp[c(13:15)]) # cols(geom,part,x,y,hole)  where geom is the MEASO area number
  r1<-r0[r0[,"part"]==1,]
  r2 <- r0[r0[,"part"]==2,] 
  r2[,"geom"]<-r2[,"geom"]+3
  r2[,"part"]<-1
  r3 <- r2  # part 2 - give this longitudes less than -180
  r3[,"x"]<-r3[,"x"]-360
  
  r<-rbind(r1,r3)
  r[,"part"]<-1
  
  sAdd180<-vect(r2,type="polygons",crs=crs(shp))
  
  sWPA<-fnJoinWParea(1,4,r,shp)
  sWPS<-fnJoinWParea(3,6,r,shp)
  sWPN<-fnJoinWParea(2,5,r,shp)
  
  dR<-values(shp[c(13:15)])
  shpNew<-vect(c(shp[c(1:12)],sWPA,sWPN,sWPS,sAdd180))
  values(shpNew)<-rbind(values(shp),dR)
  
  vR1<-shpNew; crs(vR1)<-""  # need crs="" to ensure can plot on my scales in geom_spatvector
  
  # demoplot(vR1$col,type="bar")
  
  g <-ggplot()
  g <- g + scale_x_continuous(limits=c(-190,170),breaks = seq(-180,180,30), labels=as.character(seq(-180,180,30))) +
    scale_y_continuous(limits=c(-83,-40),breaks = seq(-80,-40,10), labels=as.character(seq(-80,-40,10))) +
    theme(panel.background =element_rect(colour="black",fill = NA)
          ,plot.background = element_rect(colour="black",fill = NA)
          ,aspect.ratio = 0.6
          ,axis.title = element_blank()
          ,legend.position = "none")
  g<-g + geom_spatvector(data=vR1,aes(fill=factor(vR1$MEASO_area,levels = unique(vR1$MEASO_area))),color="black") #mapping=aes(fill=vR1$col),
  g<-g + scale_fill_manual(values = MEASOcols)
  
  return(g)
}

plotLegend<-function(dIn,oD,oL,oT,oM,oS,Xlim,Ylim,Colours){ # based on plotMEASOareaSummary
  
  BottomTickLabelOffset<- 9  # add to bottom of Ylim
  TopTickLabelOffset   <- 4  # add to top of Ylim
  LeftOffset           <-10

  Yrange<-(Ylim[2]-Ylim[1])
  
  # breakup input data ####
  Depth  <- dIn$Depth
  Light  <- dIn$Light
  Temp   <- dIn$Temp
  MLD    <- dIn$MLD
  SeaIce <- dIn$SeaIce
  
  # reduce icon width size by multiplier ####
  # icon width multiplier
  im <- 0.7
  areaXmaxAdj<-0.05
  Xlim[2]<-Xlim[2]/im  # spread out x-axis but keep icons same

  oD$rect$xmax     <- oD$rect$xmax-areaXmaxAdj
  oD$XpropPerUnit  <- oD$XpropPerUnit*im
  oL$sun$r0        <- oL$sun$r0*im
  oL$rect$w        <- oL$rect$w*im
  oT$bulb$r0       <- oT$bulb$r0*im
  oT$rect$w        <- oT$rect$w*im
  oM$wave$w        <- oM$wave$w*im
  oM$rect$w        <- oM$rect$w*im
  oS$rect$xmax     <- oS$rect$xmax-areaXmaxAdj
  oS$XpropPerUnit  <- oS$XpropPerUnit*im
  
  # depth polygon calculations ####
  Xrange<- 10 # fixed for legend
  oD$box<-list(xmin = (oD$rect$xmax-Xrange*oD$XpropPerUnit)*Xlim[2]
               ,xmax = oD$rect$xmax*Xlim[2]
               ,ymin = oD$rect$ymin*Ylim[2]
               ,ymax = (oD$rect$ymin+oD$rect$h)*Ylim[2])
  
  useProp<-Depth$prop>0
  X0    <- oD$Xstart
  
  dX<-c(X0,Depth$prop[useProp]*Depth$TotalArea)
  OriginDepth<-if(Depth$MaxHeight>=0) 0 else Depth$MaxHeight 
  dY<-c(OriginDepth,Depth$depth[useProp])   
  
  dX<-c(dX,dX[1]); dY<-c(dY,dY[length(dY)]) # returning the line to close the polygon
  dX <-  (dX-oD$Xstart)/oD$unit*oD$XpropPerUnit*Xlim[2]+oD$box$xmin
  
  
  dY<-(oD$rect$ymin+oD$rect$h)*Ylim[2]-dY/Depth$depth[length(Depth$depth)]*oD$rect$h*Yrange
  
  xTickLast<-  floor((Depth$TotalArea-oD$Xstart)/oD$unit)
  xTickEnd <- (Depth$TotalArea--oD$Xstart)/oD$unit
  
  DepthInt<-seq(0,Depth$depth[length(Depth$depth)],-1000)
  oD$box<-c(oD$box,list(
    ticksX =  c(seq(0,xTickLast,1),xTickEnd)*oD$XpropPerUnit*Xlim[2]+oD$box$xmin 
    ,ticksY = oD$box$ymax-DepthInt/(Depth$depth[length(Depth$depth)])*oD$rect$h*Yrange
    ,poly = data.frame(x= dX, y= dY)
  ))
  
  # light calculations ####
  oL$sun$x<-oL$sun$x0*Xlim[2]
  oL$sun$y<-oL$sun$y0*Ylim[2]
  oL$sun$r<-oL$sun$r0*Xlim[2]
  
  thetaInt<-2*pi/oL$sun$rayN
  theta<-seq(0,2*pi,thetaInt)
  
  oL$sun$poly<-as.data.frame(matrix(unlist(sapply(theta,function(t,sun,L,tI){
    r<-c((sun$r-L),sun$r)
    theta<-c(t,t+tI/2)
    x<-cos(theta)*r+sun$x
    y<-sin(theta)*r+sun$y
    return(rbind(theta,x,y))
  },oL$sun,oL$sun$rayLen*Xlim[2],thetaInt)),ncol=3,byrow=TRUE))
  names(oL$sun$poly)<-c("theta","x","y")
  
  # background box to hide rays
  oL$rectBG$xmin <- oL$sun$x-oL$rect$w*Xlim[2]/2
  oL$rectBG$xmax <- oL$sun$x+oL$rect$w*Xlim[2]/2
  oL$rectBG$ymin <- oL$rect$ymax-oL$sun$rayLen*Xlim[2] - oL$rect$h*Yrange
  
  # data box
  oL$rect$xmin <- oL$sun$x-oL$rect$w*Xlim[2]/2
  oL$rect$xmax <- oL$sun$x+oL$rect$w*Xlim[2]/2
  oL$rect$ymax <- oL$sun$y - (oL$sun$r^2-(oL$rect$w*Xlim[2]/2)^2)^0.5
  oL$rect$ymin <- oL$rect$ymax - oL$rect$h*Yrange
  oL$rect$ymaxBG <- oL$sun$y - ((oL$sun$r-oL$sun$rayLen*Xlim[2])^2-(oL$rect$w*Xlim[2]/2)^2)^0.5
  
  
  oL$ticksY <- oL$rect$ymax-(seq(oL$Lmin,oL$Lmax,length.out=(oL$Ints+1))-oL$Lmin)/(oL$Lmax-oL$Lmin)*oL$rect$h*Ylim[2]
  oL$ticksX <- oL$rect$xmin- c(0,2)#oL$rect$w*Xlim[2]*c(0.4,0.6)
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
  
  oT$ticksY <- (seq(oT$Tmin,oT$Tmax,length.out=(oT$Ints+1))-oT$Tmin)/(oT$Tmax-oT$Tmin)*oT$rect$h*Ylim[2]+oT$rect$ymin
  oT$ticksX <- oT$rect$xmin- c(0,2)

  oT$tempShape <-data.frame(x = c(oT$rect$xmin,oT$rect$xmin,oT$rect$xmax,oT$rect$xmax)
                            ,y = c(oT$rect$ymin
                                   ,(Temp$early-oT$Tmin)/(oT$Tmax-oT$Tmin)*oT$rect$h*Ylim[2]+oT$rect$ymin
                                   ,(Temp$late-oT$Tmin)/(oT$Tmax-oT$Tmin)*oT$rect$h*Ylim[2]+oT$rect$ymin
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
    
    boxX<-c(polyX,pXe,pX1) # polyX is a long vector
    boxY<-c(polyY,rectYbottom,rectYbottom)
    
    list(box = data.frame(x=boxX*Xlim[2],y=boxY*Ylim[2])
         ,boxXcentre = (pXe+pX1)/2*Xlim[2]
         ,polybase = data.frame(x=c(polyX,pX1),y=c(polyY,rectYtop))
         ,poly = data.frame(x=c(pX1,polyX,pXe)*Xlim[2],y=c(pY1,polyY,pYe)*Ylim[2])
         ,bottom=list(x=c(pX1,pXe)*Xlim[2],y=c(rectYbottom,rectYbottom)*Ylim[2])
         ,ticksX = pX1*Xlim[2]- c(0,2)
         ,ticksY = (rectYtop-((seq(oM$MLmin,oM$MLmax,length.out=(oM$Ints+1))-oM$MLmin)/(oM$MLmax-oM$MLmin)*oM$rect$h))*Ylim[2]
    )
  })

    # Sea ice concentration and extent (plotted over depth) ####
  sXrange<- if(is.na(oS$unit)) (log10(SeaIce$MEASOarea)-oS$Xstart) else (SeaIce$MEASOarea-oS$Xstart)/oS$unit
  
  oS$box<-list( xmin = (oS$rect$xmax-sXrange*oS$XpropPerUnit)*Xlim[2]
                ,xmax = oS$rect$xmax*Xlim[2]
                ,ymin = oS$rect$ymin*Yrange+Ylim[1]
                ,ymax = (oS$rect$ymin+oS$rect$h)*Yrange+Ylim[1])
  
  s95<-list(early = (SeaIce$conc95$early/100*SeaIce$MEASOarea-oS$Xstart)/oS$unit
            ,late = (SeaIce$conc95$late/100*SeaIce$MEASOarea-oS$Xstart)/oS$unit)
  s75<-list(early = (SeaIce$conc75$early/100*SeaIce$MEASOarea-oS$Xstart)/oS$unit
            ,late = (SeaIce$conc75$late/100*SeaIce$MEASOarea-oS$Xstart)/oS$unit)
  s25<-list(early = (SeaIce$conc25$early/100*SeaIce$MEASOarea-oS$Xstart)/oS$unit
            ,late = (SeaIce$conc25$late/100*SeaIce$MEASOarea-oS$Xstart)/oS$unit)
  oS$box<-c(oS$box,list(
    poly95 = data.frame(x=c(oS$box$xmin,oS$box$xmin
                            ,oS$box$xmin+s95$early*oS$XpropPerUnit*Xlim[2]
                            ,oS$box$xmin+s95$late*oS$XpropPerUnit*Xlim[2])
                        ,y=c(oS$box$ymin,oS$box$ymax,oS$box$ymax,oS$box$ymin))
    ,poly75 = data.frame(x=c(oS$box$xmin,oS$box$xmin
                             ,oS$box$xmin+s75$early*oS$XpropPerUnit*Xlim[2]
                             ,oS$box$xmin+s75$late*oS$XpropPerUnit*Xlim[2])
                         ,y=c(oS$box$ymin,oS$box$ymax,oS$box$ymax,oS$box$ymin))
    ,poly25 = data.frame(x=c(oS$box$xmin,oS$box$xmin
                             ,oS$box$xmin+s25$early*oS$XpropPerUnit*Xlim[2]
                             ,oS$box$xmin+s25$late*oS$XpropPerUnit*Xlim[2])
                         ,y=c(oS$box$ymin,oS$box$ymax,oS$box$ymax,oS$box$ymin))
  ))
  
  # Plot ####
  Ylim[1]<-Ylim[1]-BottomTickLabelOffset # correct Ylim after all the calcs
  Ylim[2]<-Ylim[2]+TopTickLabelOffset
  Xlim[1]<-Xlim[1]-LeftOffset
  
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
  
  AxisTickLabelSize<-4
  AxisLabelSize<-5
  # depth ####
  p<-p+geom_rect(aes(xmin=oD$box$xmin,xmax=oD$box$xmax,ymin=oD$box$ymin,ymax=oD$box$ymax),fill=Colours$depth$bg,linetype=1,colour="black",show.legend=FALSE)
  p<-p+geom_polygon(data=oD$box$poly,aes(x=x,y=y),fill=Colours$depth$main, show.legend=FALSE)
  for(i in c(1:length(oD$box$ticksX))) p<-p+geom_line(aes(x=x,y=y),data=data.frame(x=c(oD$box$ticksX[i],oD$box$ticksX[i]),y=c(oD$box$ymin-2,oD$box$ymin)),linetype=1,colour="black")
  for(i in c(1:length(oD$box$ticksY))) p<-p+geom_line(aes(x=x,y=y),data=data.frame(x=c(oD$box$xmin-2,oD$box$xmin), y=c(oD$box$ticksY[i],oD$box$ticksY[i])),linetype=1,colour="black")
  p<-p+ annotate(geom="text", x= oD$box$xmin-3, y=oD$box$ymin , label="-4000",
                 color="black", hjust=1,vjust=0.5,size=AxisTickLabelSize)
  p<-p+ annotate(geom="text", x= oD$box$xmin-3, y=oD$box$ymax , label="0",
                 color="black", hjust=1,vjust=1,size=AxisTickLabelSize)
  p<-p+ annotate(geom="text", x= oD$box$xmin-23, y=(oD$box$ymin+oD$box$ymax)/2+0.5 , label="Seafloor",
                 color="black", hjust=0,vjust=0,size=AxisLabelSize,angle=00)
  p<-p+ annotate(geom="text", x= oD$box$xmin-23, y=(oD$box$ymin+oD$box$ymax)/2-0.5 , label="Depth (m)",
                 color="black", hjust=0,vjust=1,size=AxisLabelSize,angle=0)
  
  p<-p+ annotate(geom="text", x= oD$box$xmin, y=oD$box$ymin-3 , label="0",
                 color="black", hjust=0.5,vjust=1,size=AxisTickLabelSize)
  p<-p+ annotate(geom="text", x= oD$box$xmax, y=oD$box$ymin-3 , label="10",
                 color="black", hjust=0.5,vjust=1,size=AxisTickLabelSize)
  p<-p+ annotate(geom="text", x= (oD$box$xmin+oD$box$xmax)/2, y=Ylim[1] , label="Area (million square km)",
                 color="black", hjust=0.5,vjust=0,size=AxisLabelSize)
  
  
  # light ####
  
  #p<-p+ geom_circle(aes(x0 = oL$sun$x, y0 = oL$sun$y, r = oL$sun$r), fill = Colours$light$main, show.legend=FALSE)
  p<-p+geom_polygon(data=oL$sun$poly,aes(x=x,y=y),fill=Colours$light$main,linetype=1,colour="black", show.legend=FALSE)
  #p<-p+geom_rect(aes(xmin=oL$rect$xmin,xmax=oL$rect$xmax,ymin=oL$rect$ymax,ymax=oL$rect$ymaxBG),fill=Colours$light$main,show.legend=FALSE)
  p<-p+geom_rect(aes(xmin=oL$rect$xmin,xmax=oL$rect$xmax,ymin=oL$rect$ymin,ymax=oL$rect$ymax),fill=Colours$light$bg,linetype=1,colour="black",show.legend=FALSE)
  p<-p+geom_polygon(data=oL$poly,aes(x=x,y=y),fill=Colours$light$main, show.legend=FALSE)
  for(i in c(1:length(oL$ticksY))) p<-p+geom_line(aes(x=x,y=y),data=data.frame(x=oL$ticksX,y=c(oL$ticksY[i],oL$ticksY[i])),linetype=1,colour="black")
  p<-p+ annotate(geom="text", x= oL$rect$xmin-10, y=0 , label="Watts.m-2.d-1",
                 color="black", hjust=0,vjust=0,size=AxisLabelSize,angle=90)
  p<-p+ annotate(geom="text", x= oL$rect$xmin-3, y=oL$rect$ymax , label=oL$Lmin,
                 color="black", hjust=1,vjust=1,size=AxisTickLabelSize)
  p<-p+ annotate(geom="text", x= oL$rect$xmin-3, y=oL$rect$ymin , label=oL$Lmax,
                 color="black", hjust=1,vjust=0,size=AxisTickLabelSize)
  p<-p+ annotate(geom="text", x= (oL$rect$xmin+oL$rect$xmax)/2, y=Ylim[1] , label="PAR",
                 color="black", hjust=0.5,vjust=0,size=AxisLabelSize)

  p<-p+geom_line(aes(x=x,y=y),data=data.frame(x=c(oL$rect$xmin,oL$rect$xmin),y=c(oL$rect$ymin,oL$rect$ymin-2)),linetype=1,colour="black")
  p<-p+geom_line(aes(x=x,y=y),data=data.frame(x=c(oL$rect$xmax,oL$rect$xmax),y=c(oL$rect$ymin,oL$rect$ymin-2)),linetype=1,colour="black")
  p<-p+ annotate(geom="text", x= oL$rect$xmin,y=oL$rect$ymin-3, label="R",
                 color="black", hjust=0.5,vjust=1,size=AxisTickLabelSize)
  p<-p+ annotate(geom="text", x= oL$rect$xmax,y=oL$rect$ymin-3, label="E",
                 color="black", hjust=0.5,vjust=1,size=AxisTickLabelSize)
  
  # thermometer ####
  p<-p+ geom_circle(aes(x0 = oT$bulb$x, y0 = oT$bulb$y, r = oT$bulb$r), fill = Colours$temp$main, show.legend=FALSE)
  p<-p+geom_rect(aes(xmin=oT$rect$xmin,xmax=oT$rect$xmax,ymin=oT$rect$ymin,ymax=oT$rect$ymax),fill=Colours$temp$bg,linetype=1,colour="black",show.legend=FALSE)
  p<-p+geom_polygon(data=oT$tempShape,aes(x=x,y=y),fill=Colours$temp$main, show.legend=FALSE)
  for(i in c(1:length(oT$ticksY))) p<-p+geom_line(aes(x=x,y=y),data=data.frame(x=oT$ticksX,y=c(oT$ticksY[i],oT$ticksY[i])),linetype=1,colour="black")

  p<-p+ annotate(geom="text", x= oT$rect$xmin-8, y=0 , label="      degrees C",
                 color="black", hjust=0,vjust=0,size=AxisLabelSize,angle=90)
  p<-p+ annotate(geom="text", x= oT$rect$xmin-3, y=oT$rect$ymax , label=oT$Tmax,
                 color="black", hjust=1,vjust=1,size=AxisTickLabelSize)
  p<-p+ annotate(geom="text", x= oT$rect$xmin-3, y=oT$rect$ymin , label=oT$Tmin,
                 color="black", hjust=1,vjust=0,size=AxisTickLabelSize)
  p<-p+ annotate(geom="text", x= (oT$rect$xmin+oT$rect$xmax)/2, y=Ylim[1] , label="Temperature   ",
                 color="black", hjust=0.5,vjust=0,size=AxisLabelSize)
  
  p<-p+geom_line(aes(x=x,y=y),data=data.frame(x=c(oT$rect$xmin,oT$rect$xmin),y=c(oT$rect$ymax,oT$rect$ymax+2)),linetype=1,colour="black")
  p<-p+geom_line(aes(x=x,y=y),data=data.frame(x=c(oT$rect$xmax,oT$rect$xmax),y=c(oT$rect$ymax,oT$rect$ymax+2)),linetype=1,colour="black")
  
  p<-p+ annotate(geom="text", x= oT$rect$xmin,y=oT$rect$ymax+3, label="R",
                 color="black", hjust=0.5,vjust=0,size=AxisTickLabelSize)
  p<-p+ annotate(geom="text", x= oT$rect$xmax,y=oT$rect$ymax+3, label="E",
                 color="black", hjust=0.5,vjust=0,size=AxisTickLabelSize)
  
  # MLD ####
  p<-p+geom_polygon(data=oM$wave$df,aes(x=x,y=y),fill=Colours$MLD$main,linetype=1,colour="black", show.legend=FALSE)
  p<-p+geom_polygon(data=oM$box,aes(x=x,y=y),fill=Colours$MLD$bg,linetype=1,colour="black", show.legend=FALSE)
  p<-p+geom_polygon(data=oM$poly,aes(x=x,y=y),fill=Colours$MLD$main, show.legend=FALSE)
  p<-p+geom_polygon(data=oM$polybase,aes(x=x,y=y),fill=Colours$MLD$main,linetype=1,colour="black", show.legend=FALSE)
  
  for(i in c(1:length(oM$ticksY))) p<-p+geom_line(aes(x=x,y=y),data=data.frame(x=oM$ticksX,y=c(oM$ticksY[i],oM$ticksY[i])),linetype=1,colour="black")
  p<-p+ annotate(geom="text", x= oM$ticksX[2]-6, y=0 , label="Depth (m)",
                 color="black", hjust=0,vjust=0,size=AxisLabelSize,angle=90)
  p<-p+ annotate(geom="text", x= oM$ticksX[2]-0.5, y=oM$ticksY[1] , label=oM$MLmin,
                 color="black", hjust=1,vjust=1,size=AxisTickLabelSize)
  p<-p+ annotate(geom="text", x= oM$ticksX[2]-0.5, y=oM$ticksY[oM$Ints+1] , label=oM$MLmax,
                 color="black", hjust=1,vjust=0,size=AxisTickLabelSize)
  p<-p+ annotate(geom="text", x= oM$boxXcentre, y=Ylim[1] , label=" Mixed Layer",
                 color="black", hjust=0.5,vjust=0,size=AxisLabelSize)

  p<-p+geom_line(aes(x=x,y=y),data=data.frame(x=c(oM$bottom$x[1],oM$bottom$x[1]),y=c(oM$bottom$y[1],oM$bottom$y[1]-2)),linetype=1,colour="black")
  p<-p+geom_line(aes(x=x,y=y),data=data.frame(x=c(oM$bottom$x[2],oM$bottom$x[2]),y=c(oM$bottom$y[1],oM$bottom$y[1]-2)),linetype=1,colour="black")
  
  p<-p+ annotate(geom="text", x= oM$bottom$x[1],y=oM$bottom$y[1]-3, label="R",
                 color="black", hjust=0.5,vjust=1,size=AxisTickLabelSize)
  p<-p+ annotate(geom="text", x= oM$bottom$x[2],y=oM$bottom$y[2]-3, label="E",
                 color="black", hjust=0.5,vjust=1,size=AxisTickLabelSize)
  
  
  # Sea ice ####
    if(sum(oS$box$poly95$x>oS$box$xmin)>0){
    # Sea ice (data are for open water)  
    # p<-p+geom_rect(aes(xmin=oS$box$xmin,xmax=oS$box$xmax,ymin=oS$box$ymin,ymax=oS$box$ymax),fill=Colours$ice$bg,linetype=1,colour="black",show.legend=FALSE)
    p<-p+geom_polygon(data=oS$box$poly95,aes(x=x,y=y),fill=Colours$ice$c95,linetype=1,colour="black", show.legend=FALSE)
    p<-p+geom_polygon(data=oS$box$poly75,aes(x=x,y=y),fill=Colours$ice$c75, show.legend=FALSE)
    p<-p+geom_polygon(data=oS$box$poly25,aes(x=x,y=y),fill=Colours$ice$c25, show.legend=FALSE)
    # y ticks
    p<-p+geom_line(aes(x=x,y=y),data=data.frame(x=c(oS$box$xmin-2,oS$box$xmin), y=c(oS$box$ymin,oS$box$ymin)),linetype=1,colour="black")
    p<-p+geom_line(aes(x=x,y=y),data=data.frame(x=c(oS$box$xmin-2,oS$box$xmin), y=c(oS$box$ymax,oS$box$ymax)),linetype=1,colour="black")
    p<-p+ annotate(geom="text", x= oS$box$xmin-3,y=oS$box$ymin, label="E",
                   color="black", hjust=1,vjust=0,size=AxisTickLabelSize)
    p<-p+ annotate(geom="text", x= oS$box$xmin-3,y=oS$box$ymax, label="R",
                   color="black", hjust=1,vjust=0.5,size=AxisTickLabelSize)
    p<-p+ annotate(geom="text", x= oS$box$xmin-20, y=(oS$box$ymin+oS$box$ymax)/2+1 , label="Spring",
                   color="black", hjust=0,vjust=0,size=AxisLabelSize,angle=0)
    p<-p+ annotate(geom="text", x= oS$box$xmin-20, y=(oS$box$ymin+oS$box$ymax)/2 , label="Sea Ice",
                   color="black", hjust=0,vjust=1,size=AxisLabelSize,angle=0)
    
    
  }
  
  
  return(p)
}

p<-plotLegend(MEASOareaSummary[[1]],oD,oL,oT,oM,oS,Xlim,Ylim,Colours)  
p

# 4.0 Test Routine ####

# plot MEASO summaries

for(i in 10:12){
p<-plotMEASOareaSummary(MEASOareaSummary[[i]],oD,oL,oT,oM,oS,Xlim,Ylim,doTicks,Colours,ContShelf=TRUE)
print(p)
}


# loading the required libraries 
library(ggplot2)
library(png) 
library(tidyterra)
library(patchwork)
library(geodata)
library(terra)


MEASOshp<-cbind(MEASOshp,as.data.frame(MEASOcols))
names(MEASOshp)<-c(names(MEASOshp)[1],"col")

baseMap<-plotMEASObaseMap(MEASOshp)
baseMap

pSize<-c(0.25,0.25)
pLoc<-list( a01 = c(0.60,0.3) # AOA
           ,a02 = c(0.60,0.72) # AON
           ,a03 = c(0.60,0.59) # AOS
           ,a04 = c(0.80,0.30) # CIA
           ,a05 = c(0.80,0.70) # CIN
           ,a06 = c(0.80,0.53) # CIS
           ,a07 = c(0.97,0.3) # EIA
           ,a08 = c(0.97,0.60) # EIN
           ,a09 = c(0.97,0.44) # EIS
           ,a10 = c(0.365,0.14) # EPA
           ,a11 = c(0.365,0.60) # EPN
           ,a12 = c(0.365,0.34) # EPS
           ,a13 = c(0.22,0.14) # WPA
           ,a14 = c(0.22,0.60) # WPN
           ,a15 = c(0.22,0.40) # WPS
)  # coordinates for bottom-right corner because size may vary

g<-baseMap
for(i in seq(1,15,1)){
  p<-plotMEASOareaSummary(MEASOareaSummary[[i]],oD,oL,oT,oM,oS,Xlim,Ylim,doTicks,Colours,ContShelf=TRUE)
  
  g<-g + inset_element( p
                ,left   = pLoc[[i]][1]-pSize[1] 
                ,bottom = pLoc[[i]][2] 
                ,right  = pLoc[[i]][1] 
                ,top    = pLoc[[i]][2]+pSize[2])
       }

#p<-plotLegend(MEASOareaSummary[[1]],oD,oL,oT,oM,oS,Xlim,Ylim,Colours)  
img <-  readPNG('MEASOlegend.png')
g<- g + annotation_raster( img
                      ,xmin   = 0.93-pSize[1]*2
                      ,xmax  = 0.93 
                      ,ymin = 0.005 
                      ,ymax    = 0.005+pSize[2])
g


