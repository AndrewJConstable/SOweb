# Net primary production summaries
# 
# Data source: ####
# Ryan-Keogh, T., Thomalla, S., Chang, N., & Moalusi, T. (2023)
# Net primary production from the Eppley-VGPM, Behrenfeld-VGPM, 
# Behrenfeld-CbPM, Westberry-CbPM and Silsbe-CAFE algorithms (1.2) [Data set]. 
# Zenodo. https://doi.org/10.5281/zenodo.10829966
# Version 1.2 18 March 2024


# This file generates a dataframe of
#  Var       : Variable
#  MEASO     : MEASO area
#  mean,sd   : total net primary productivity for year (1 July to 30 June) 

# 1. Notes ####
#   * Regular grid with longitude range -180:180 (Var: rsds)
#   * Cell centres are inside range with cell edges at range limits
#   * MEASO shapes are projected with longitude range -180-180
#   * Each ACCESS cell may contribute to more than one MEASO cell - proportional contribution to MEASO calculated
#   * units of NPP are "mg C m^-2 d^-1"  - thus multiply NPP x cell area  x days

library(RNetCDF)
library(terra)
library(viridis)
library(lubridate)
library(ggplot2)
library(patchwork)

# 2. Input Data ####
DoInFull<-FALSE

#         2.1 Bounds of ACCESS data grid to use ####
domainBounds<-list(limLon=c(-180,180),limLat=c(-80,-30),resolution=0.1) # also for regular grid to use
rExt<-ext(c(domainBounds$limLon,domainBounds$limLat))
#         2.2 Directories ####
rootData<-"/Users/acon/Desktop/_w/_d/Ryan-Keogh 2023 NPP/"
rootMEASOshp<-"/Users/acon/Desktop/_w/_d/Shapefiles/MEASO polygons/"

#         2.3 Files ####

fMEASOshp<-"MEASO_polygons.shp"

fGen<-list(root   = "/Users/acon/Desktop/_w/_d/Ryan-Keogh 2023 NPP/"
                 ,ext    = ".nc"
                 ,miss   = NaN
                 ,res    = c(0.25,0.25)
                ) # end detail
fVar<-list(B_CbPM   = list(fname = "1998_2023_NPP_BEHRENFELD_CBPM_25KM_8D" 
                        ,var   = "Behrenfeld-CbPM"  #lon, lat, time
                        ) # end B_CbPM
          ,B_VGPM   = list(fname = "1998_2023_NPP_BEHRENFELD_VGPM_25KM_8D"
                        ,var   = ""
                        ) # end B_VGPM
          ,E_VGPM   = list(fname = "1998_2023_NPP_EPPLEY_VGPM_25KM_8D"
                           ,var   = ""
                          ) # end E_VGPM
          ,S_CAFE   = list(fname = "1998_2023_NPP_SILSBE_CAFE_25KM_8D"
                           ,var   = ""
                           ) # end S_CAFE
          ,W_CbPM   = list(fname = "1998_2023_NPP_WESTBERRY_CBPM_25KM_8D"
                           ,var   = ""
                           ) # end W_CbPM
) # end fVar

# exploratory code ####
ignore<-TRUE
if(!ignore){
ncin<-open.nc(paste0(rootData,fNPP_B_CBPM,".nc"))
lat<-var.get.nc(ncin,"lat") 
lon<-var.get.nc(ncin,"lon")
t<-var.get.nc(ncin,"time")
}#####

Day0<-as.Date("1998-01-01")

# extract data for following
eYears   <- c(2001:2020)
ePeriods <- list(Summer = c(1,3)
                ,Autumn = c(4,6)
                ,Winter = c(7,9)
                ,Spring = c(10,12))

# cropping original datasets to Southern Ocean ####

if(DoInFull){
lapply(fVar,function(f,g,e){
  print(f$fname)
  r<-rast(paste0(g$root,f$fname,".nc"))
  r<-crop(r,e)
  writeRaster(r, paste0(g$root,f$fname,"_MEASO.tif"))
},fGen,rExt)
}

# load all raster bricks and reduce to quarterly means

rAll<-lapply(fVar,function(f,g){rast(paste0(g$root,f$fname,"_MEASO.tif"))},fGen)

tSteps<-c(1:length(rAll[[1]][1]))
Dates  <- Day0+(tSteps*8-8)
Years  <- year(Dates)
Months <- month(Dates)

rAll<-lapply(rAll,function(r,y,m,ey,ep){
            rast(lapply(ey,function(ey,r,y,m,ep){
              rast(lapply(ep,function(p,r,m){
                  mean(r[[m%in%p]],na.rm=TRUE)
                  },r[[y==ey]],m[y==ey]))
            },r,y,m,ep))
            },Years,Months,eYears,ePeriods)

# area of each cell

cs<-cellSize(rAll[[1]][[1]],unit="km")

# create dataframe rows = rasters; cols 1= year, 2= number of days in each quarter

QrtrDays<-do.call(rbind,lapply(eYears,function(ey,y,m,ep){
  do.call(rbind,lapply(ep,function(p,Year,m){
    data.frame(Year,Days=sum(m%in%p)*8) # calc number of days
  },ey,m[y==ey]))
},Years,Months,ePeriods))
QrtrDays<-data.frame(Year=QrtrDays$Year,Season=rep(names(ePeriods),length(eYears)),Days=QrtrDays$Days)

# generate cell means for each time period * cell area (km2) x 1E-6  units become C 1E3 t d-1 in cell

rMn<-rast(lapply(c(1:length(rAll[[1]][1])),function(t,r,cs,qd){
          mean(rast(lapply(r,function(r,t){r[[t]]},t)))*cs*qd[t]*1E-6
          },rAll,cs,QrtrDays$Days))


ignore<-TRUE
if(!ignore){
  # exploratory code
r<-rast(paste0(fGen$root,fVar[[1]]$fname,"_MEASO.tif"))
plot(r)
plot(mshp)
extract(r)
}

mshp<-vect(paste0(rootMEASOshp,fMEASOshp))
NPPyr<-do.call(rbind,lapply(eYears,function(y,r,Y,mshp){
  ry<-sum(r[[Y==y]],na.rm=TRUE)
  df<-extract(ry,mshp,weight=TRUE)
  res<-aggregate((df$sum*df$weight),by=list(MEASO=df$ID),FUN=sum,na.rm=TRUE)
  return(data.frame(MEASO=res$MEASO,Year=rep(y,nrow(res)),NPP=res$x))
},rMn,QrtrDays$Year,mshp))

NPPmeaso<-merge(aggregate(NPPyr$NPP,by=list(MEASO=NPPyr$MEASO),FUN=mean)
           ,aggregate(NPPyr$NPP,by=list(MEASO=NPPyr$MEASO),FUN=sd)
           ,by="MEASO")

#                1      2      3      4      5      6      7      8      9      10    11     12     13     14    15
MEASOnames<-c("AOA", "AON", "AOS", "CIA", "CIN", "CIS", "EIA", "EIN", "EIS", "EPA", "EPN", "EPS", "WPA", "WPN", "WPS")

NPPmeaso<-cbind(MEASOnames,NPPmeaso)
names(NPPmeaso)<-c("MEASO","ID","Mean","SD")
NPPmeaso$MEASO<-factor(NPPmeaso$MEASO,levels=c("AOA", "AOS", "AON"
                                     , "CIA", "CIS", "CIN"
                                     , "EIA", "EIS", "EIN" 
                                     , "EPA", "EPS", "EPN"
                                     , "WPA", "WPS", "WPN"))

p<-ggplot(NPPmeaso,aes(MEASO,Mean,colour=MEASO))+geom_bar(stat="identity")+geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD,colour=MEASO))
p



Sectors<-c("AO","CI","EI","WP","EP")
useSA<-do.call(rbind,lapply(c(1:5),function(i,s,n){
  r1<-which(substr(n,1,2)%in%s[i] & substr(n,3,3)!="N")
  return(cbind(rep(i,length(r1)),r1))},Sectors,MEASOnames))
useSectorAreas<-rep(0,15)
useSectorAreas[useSA[,2]]<-useSA[,1]

NPPyr<-cbind(NPPyr,Sector=rep(useSectorAreas,20))
NPPyrSector<-aggregate(NPPyr$NPP,by=list(Year=NPPyr$Year,Sector=NPPyr$Sector),FUN=sum)
NPPyrSector<-NPPyrSector[NPPyrSector$Sector!=0,]
names(NPPyrSector)<-c("Year","Sector","NPP")
NPPsector<-merge(aggregate(NPPyrSector$NPP,by=list(Sector=NPPyrSector$Sector),FUN=mean)
                ,aggregate(NPPyrSector$NPP,by=list(Sector=NPPyrSector$Sector),FUN=sd)
                ,by="Sector")

NPPsector<-cbind(Sectors,NPPsector)
names(NPPsector)<-c("Sector","ID","Mean","SD")
NPPsector$Sector<-factor(NPPsector$Sector,levels=Sectors)

p<-ggplot(NPPsector,aes(Sector,Mean,colour=Sector))+geom_bar(stat="identity")+geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD,colour=Sector))
p


NPPtotal<-aggregate(NPPyr$NPP,by=list(Year=NPPyr$Year),FUN=sum)
print(paste0("Total NPP: mean = ",mean(NPPtotal$x),", SD = ",sd(NPPtotal$x)))

#################################
# Bubble plots on maps

# Use basemap from MEASO map plot variables.R

# bubble locations
dNPPa<-data.frame(
  Group = rep("Area",15)
, MEASO = c("AOA","AON","AOS","CIA","CIN","CIS","EIA","EIN","EIS","EPA","EPN","EPS","WPA","WPN","WPS")
, x =     c( 0.48, 0.48, 0.48, 0.72, 0.72, 0.72, 0.95, 0.95, 0.95, 0.22, 0.22, 0.22,  0.05, 0.05, 0.05)
, y =     c( 0.40, 0.93, 0.75, 0.40, 0.93, 0.70, 0.40, 0.80, 0.60, 0.25, 0.80, 0.50,  0.3, 0.80, 0.60)
)
dNPPa<-merge(dNPPa,NPPmeaso,by="MEASO")

dNPPs<-data.frame(
    Group = rep("Sectors minus Northern",5)
  , MEASO = c("AO", "CI", "EI", "EP", "WP")
  , x =     c(0.48, 0.72, 0.95, 0.22, 0.05)
  , y =     c(0.08, 0.08, 0.08, 0.08, 0.08)
)
dNPPs<-merge(dNPPs,NPPsector,by.x="MEASO",by.y="Sector")

dNPP<-rbind(dNPPa,dNPPs)
dNPP$Group<-factor(dNPP$Group)
dNPP<-cbind(dNPP,NPP=dNPP$Mean*1E-6,CV=dNPP$SD/dNPP$Mean)
# do bubble plot
p<-ggplot(dNPP,aes(x=x,y=y,size=NPP, fill=Group))
p<-p+ geom_point(alpha=1,shape=21, color="black")  # alpha is opacity
p<-p+geom_text(aes(label=MEASO),size=2.8)
p<-p+ scale_size_area(name="NPP"  #,guide="legend" ,position="bottom" #,range  = c(0, 0.7)
                      ,max_size=30, breaks = c(0.1, 0.3, 0.6))
p<-p+  scale_fill_manual(values=c("white","grey"),aesthetics="fill")
p<-p+ coord_fixed(ratio=0.6) + xlim(0,1) + ylim(0,1) + xlab("")+ylab("")+
  theme(
    axis.text.x = element_blank()
    ,axis.text.y = element_blank()
    ,axis.ticks = element_blank()
    ,panel.grid.major = element_blank()
    , panel.grid.minor = element_blank()
    ,panel.background = element_rect(fill='transparent')
    ,plot.background = element_rect(fill='transparent', color=NA) #transparent plot bg
    ,legend.position = "bottom"
   # ,legend.position.inside = c(0.8, 0.1)
    ,legend.background = element_rect(fill='transparent', color=NA)
     ) # end theme
p<-p+guides(size=guide_legend(position="bottom",title="NPP million tonnes")
          , fill=guide_legend(position="top",title="Group",override.aes = list(size = 5),theme =theme(element_text(size=20))))
p

g<-baseMap
g<-g+inset_element( p
                    ,left   = -0.1 
                    ,bottom = -0.3
                    ,right  = 1.1 
                    ,top    = 1.1
                    ,on_top = TRUE)
g



#plot(vCellsO, y="area")

tmp<-vCellsO
tmpWts<-wtsO

Ncols<-50
tmpRamp<-viridis(Ncols+1)
valWts<-as.vector(tmpWts)
valWts<-valWts[!(valWts==0 | is.na(valWts))]
AreaMin<-min(valWts)
AreaMax<-max(valWts)

plot(MEASOshp)
tmp$area[tmp$area==0]<-NA
tmpCols<-tmpRamp[round((tmp$area-AreaMin)/(AreaMax-AreaMin)*Ncols+1)]
polys(tmp, col=tmpCols, border=NA)

plot(MEASOshp)
for(i in c(1:15)) {
  tmp$area<-tmpWts[,i]
  tmp$area[tmp$area==0]<-NA
#  plot(tmp,"area")
  tmpCols<-tmpRamp[round((tmp$area-AreaMin)/(AreaMax-AreaMin)*Ncols+1)]
#  plot(MEASOshp,main=paste0("MEASO ",MEASOnames[i]))
  polys(tmp, col=tmpCols, border=NA)
}
