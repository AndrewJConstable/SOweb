# ACCESS ESM1.5 SSP5-8.5
# supported by Tilo Ziehn

# Variable definitions:
# https://clipc-services.ceda.ac.uk/dreq/mipVars.html

# extractions undertaken using code in MEASO_ACCESS_ESM_15_SSP585.R

# This file generates a dataframe of
# 
#  MEASO     : MEASO area
#  Depth     : LT2000m, DT2000m, All
#  Var1
#  Var2
#  ...
#  Var.n

# 1. Notes ####
#   * ACCESS ocean and ice models use irregular grid with longitude range 0-360 (Vars: thetao, siconc, mlotst)
#   * ACCESS atmosphere model uses regular grid with longitude range 0-360 (Var: rsds)
#   * Cell centres are on integer degrees resulting in cells crossing the boundaries at 0 and at 180/-180
#   * MEASO shapes are projected with longitude range -180-180
#   * Each ACCESS cell may contribute to more than one MEASO cell.  
#   * The weight of contribution of an ACCESS cell value to a MEASO area statistic is 
#         the proportion of the ACCESS cell in the MEASO area multiplied by the area of the ACCESS cell
#   * This routine generates weights for each ACCESS cell in a MEASO area 
#          to be used for generating statistics in each MEASO area
#          i.e. for each MEASO area, a vector of weights corresponding to the ACCESS grid
#   * The weights are saved in a dataframe - rows(ACCESS cells), cols(MEASO areas)

library(RNetCDF)
library(terra)
library(viridis)
library(vegan)
library(ggplot2)

# 2. Input Data ####

#         2.1 Bounds of ACCESS data grid to use ####
domainBounds<-list(limLon=c(0,360),limLat=c(-80,-30),resolution=0.1) # also for regular grid to use

#         2.2 Directories ####
rootData<-"/Users/acon/Desktop/_w/_d/ACCESS-ESM1.5\ SSP585/"
rootMEASOshp<-"/Users/acon/Desktop/_w/_d/Shapefiles/MEASO polygons/"
rootTimeSeries<-"/Users/acon/Desktop/_w/_d/ACCESS-ESM1.5\ SSP585/timeseries/"

#         2.3 Files ####
fCellAreaO     <- "areacello_Ofx_ACCESS-ESM1-5_ssp585_r10i1p1f1_gn.nc" # Cell Areas - ocean/ice grid
fCellAreaAtmos <- "areacella_fx_ACCESS-ESM1-5_ssp585_r10i1p1f1_gn.nc"

fMEASOshp<-"MEASO_polygons.shp"

fVar<-list(root   = "/Users/acon/Desktop/_w/_d/ACCESS-ESM1.5\ SSP585/timeseries/"
          ,ext    = ".Rdata"
          ,temp   = list(var   = "thetao"
                        ,y1900 = "1900-1919"
                        ,y1920 = "1920-1939"
                        ,y1940 = "1940-1959"
                        ,y1960 = "1960-1979"
                        ,y1980 = "1980-1999"
                        ,y2000 = "2000-2019"
                        ,y2020 = "2020-2039"
                        ,y2040 = "2040-2059"
                        ,y2060 = "2060-2079"
                        ,y2080 = "2080-2099")
          ,SeaIce = list(var   = "siconc"
                         ,y1900 = "1900-1919"
                         ,y1920 = "1920-1939"
                         ,y1940 = "1940-1959"
                         ,y1960 = "1960-1979"
                         ,y1980 = "1980-1999"
                         ,y2000 = "2000-2019"
                         ,y2020 = "2020-2039"
                         ,y2040 = "2040-2059"
                         ,y2060 = "2060-2079"
                         ,y2080 = "2080-2099")
          ,MLD    = list(var   = "mlotst"
                         ,y1900 = "1900-1919"
                         ,y1920 = "1920-1939"
                         ,y1940 = "1940-1959"
                         ,y1960 = "1960-1979"
                         ,y1980 = "1980-1999"
                         ,y2000 = "2000-2019"
                         ,y2020 = "2020-2039"
                         ,y2040 = "2040-2059"
                         ,y2060 = "2060-2079"
                         ,y2080 = "2080-2099")
          ,light  = list(var   = "rsds"
                         ,y1900 = "1900-1919"
                         ,y1920 = "1920-1939"
                         ,y1940 = "1940-1959"
                         ,y1960 = "1960-1979"
                         ,y1980 = "1980-1999"
                         ,y2000 = "2000-2019"
                         ,y2020 = "2020-2039"
                         ,y2040 = "2040-2059"
                         ,y2060 = "2060-2079"
                         ,y2080 = "2080-2099")
) # end fVar

#         2.4 Variable and month Data files to use ####

useVar<-c("light","MLD","temp","SeaIce")
useMonth <-c(11,12)
DepthBreak<-2000  # NA means return only ALL

# 3. Functions ####
#         3.1 weighted percentiles from vectors of data, weights ####
PCweighted<-function(pc,d,wt){  # vectors - percentiles, data, weights + whether to return minimum and maximum as well
  sdf      <-data.frame(d = d,wt = wt)
  delNrows <- sum(is.na(sdf[,"d"]) | sdf[,"wt"]==0)
  if(delNrows==nrow(sdf)) return(rep(NA,length(pc)))
  if(delNrows==(nrow(sdf)-1)) return(rep(sdf[!(is.na(sdf[,"d"]) | sdf[,"wt"]==0),"d"],length(pc)))
  
  sdf      <- sdf[!(is.na(sdf[,"d"]) | sdf[,"wt"]==0),]
  sdf      <- sdf[sdf[,"wt"]>0,]
  sdf      <- sdf[order(sdf$d),] # sort by data in ascending order

  cumWt    <- sapply(seq(1,nrow(sdf),1),function(i,df){sum(df[c(1:i),"wt"])},sdf)
  propWt <- cumWt/sum(sdf$wt) # calculate cumulative area and make proportion
  sdf      <- cbind(sdf,propWt)
  
  #      read area percentiles according to vector of percentiles
  #               (determine index of row given each percentile )
  #      if pc falls between two values, the higher value is returned
  #      eliminate zero weights
  
  res1<-sapply(seq(1,length(pc),1),function(i,pc,df){
    return(df[(sum(pc[i]>df$propWt)+1),"d"])
  },pc,sdf)
  names(res1)<-NULL
  return(res1)
} # end function

MeanWeighted<-function(d,wt){  # vectors - percentiles, data, weights + whether to return minimum and maximum as well
  sdf      <-data.frame(d = d,wt = wt)
  delNrows <- sum(is.na(sdf[,"d"]) | sdf[,"wt"]==0)
  if(delNrows==nrow(sdf)) return(NA)
  if(delNrows==(nrow(sdf)-1)) return(sdf[!(is.na(sdf[,"d"]) | sdf[,"wt"]==0),"d"])
  sdf      <- sdf[!(is.na(sdf[,"d"]) | sdf[,"wt"]==0),]
  res1<-weighted.mean(sdf[,"d"], sdf[,"wt"])
  names(res1)<-NULL
  return(res1)
} # end function


#         3.2 generate cell weights for ACCESS grids ####
# returns dataframe: rows(ACCESS cells - lon-lat coords are according to dataframe of ACCESS grid)
#                    cols(MEASO areas)
fnWtsOcean<-function(){

# extract details then mask using domainBounds ####
ncinO<-open.nc(paste0(rootData,fCellAreaO))
latO<-var.get.nc(ncinO,"latitude",start=c(NA,NA),count=c(NA,NA)) # note 2-D array cols(lat) rows(lon) lat[,1]=-77.8766251, lat[,300]=65.2104416
lonO<-var.get.nc(ncinO,"longitude",start=c(NA,NA),count=c(NA,NA)) # note 2-D array cols(lat) rows(lon) - also lon[1,]=80.5 lon[360,]=79.5
maskLonO<- lonO[,1]>=domainBounds$limLon[1] & lonO[,1]<domainBounds$limLon[2]
maskLatO<- latO[1,]>=domainBounds$limLat[1] & latO[1,]<domainBounds$limLat[2]
nLat<-sum(maskLatO)
nLon<-sum(maskLonO)

dlonO<-lonO[maskLonO,maskLatO]
dlatO<-latO[maskLonO,maskLatO]

CellAreaO<-var.get.nc(ncinO,"areacello",start=c(NA,NA),count=c(NA,NA)) # matrix rows(longitude) cols (latitude)
CellO_bnds_lon<-var.get.nc(ncinO,"vertices_longitude",start=c(NA,NA,NA),count=c(NA,NA,NA))  # matrix rows(bounds) cols (longitude) - start 0 end 360
CellO_bnds_lat<-var.get.nc(ncinO,"vertices_latitude",start=c(NA,NA,NA),count=c(NA,NA,NA))  # matrix rows(bounds) cols (latitude) - start minLat end maxLat

CellAreaO<-CellAreaO[maskLonO,maskLatO]
xminO<-CellO_bnds_lon[1,maskLonO,maskLatO]
xmaxO<-CellO_bnds_lon[2,maskLonO,maskLatO]
yminO<-CellO_bnds_lat[2,maskLonO,maskLatO]
ymaxO<-CellO_bnds_lat[3,maskLonO,maskLatO]

# dataframe from ACCESS ocean grid - vector of weights should correspond with rows in this grid

dfO<-data.frame(lon=as.vector(dlonO),lat=as.vector(dlatO),area=as.vector(CellAreaO),xmin=as.vector(xminO),xmax=as.vector(xmaxO),ymin=as.vector(yminO),ymax=as.vector(ymaxO))

#################################################
# change longitude to -180:180.
# current range for longitude is 0.5:359.5  with 1 degree resolution
# note that longitude in cell bounds has xmin = (0:359) & xmax = (1:0)
#  the maximum bound should be 360.

tMask<-dfO[,"lon"]>180
dfO[tMask,"lon"]<-(dfO[tMask,"lon"]-360)
dfO[dfO[,"xmin"]>=180,"xmin"]<-dfO[dfO[,"xmin"]>=180,"xmin"]-360
dfO[tMask & dfO[,"xmax"]>=180,"xmax"]<-dfO[tMask & dfO[,"xmax"]>=180,"xmax"]-360
dfO[dfO[,"xmax"]==(-360),"xmax"]<-0

#if(sum(tMask & dfO[,"xmin"]>=180)>0) dfO[tMask & dfO[,"xmin"]>=180,"xmin"]<-dfO[tMask & dfO[,"xmin"]>=180,"xmin"]-360
#if(sum(tMask & dfO[,"xmax"]>=180)>0) dfO[tMask & dfO[,"xmax"]>=180,"xmax"]<-dfO[tMask & dfO[,"xmax"]>=180,"xmax"]-360


# create SpatVector of polygons
vCellsO<-vect(lapply(c(1:nrow(dfO)),function(i,df){
       p<-cbind(c(df[i,"xmin"],df[i,"xmin"],df[i,"xmax"],df[i,"xmax"],df[i,"xmin"])
               ,c(df[i,"ymin"],df[i,"ymax"],df[i,"ymax"],df[i,"ymin"],df[i,"ymin"]))
      return(vect(p,type="polygons",crs="epsg:4326"))
       }
       ,dfO))
values(vCellsO)<-data.frame(area=dfO[,"area"],expanse=expanse(vCellsO),cell=seq(1,length(vCellsO),1))
# plot(vCellsO)

MEASOshp<-vect(paste0(rootMEASOshp,fMEASOshp))

# use to check if needed
# isRelated<-is.related(vCellsO,MEASOshp[1], "intersects")
# plot(vCellsO[isRelated])


MEASOcellWtsO<-do.call(cbind,lapply(c(1:length(MEASOshp)),function(i,m,vC){
  isRelated<-is.related(vC,m[i], "intersects")
  ret<-rep(0,length(isRelated))
  iC<-vC[m[i]]
# plot(iC,y=iC$area)
  iC<-intersect(m[i],vC)
  aC<-expanse(iC)/iC$expanse*iC$area
  ret[iC$cell]<-aC  #[order(iC$cell)]
  return(ret)
  },MEASOshp,vCellsO))

check<-FALSE
if(check){
  tmp<-data.frame(lon=dfO[,"lon"],lat=dfO[,"lat"],cell=MEASOcellWtsO[,1])
  r<-vect(tmp[tmp[,"cell"]>0,])
  plot(r,col=ceiling(r$cell*1E-9))
} # end if
return(MEASOcellWtsO)
} # end function - ocean grid 

fnWtsAtmos<-function(){
  # extract details then mask using domainBounds
  ncin<-open.nc(paste0(rootData,fCellAreaAtmos))
  lat<-var.get.nc(ncin,"lat",start=c(NA),count=c(NA))
  lon<-var.get.nc(ncin,"lon",start=c(NA),count=c(NA))
  maskLat<-lat>=domainBounds$limLat[1] & lat<=domainBounds$limLat[2]
  maskLon<-lon>=domainBounds$limLon[1] & lon<=domainBounds$limLon[2]
  
  nLat<-length(lat[maskLat])
  nLon<-length(lon[maskLon])
  
  dlat<-matrix(lat[maskLat],nrow=nLon,ncol=nLat,byrow=TRUE)
  dlon<-matrix(lon[maskLon],nrow=nLon,ncol=nLat)
  dlon[dlon>180]<-(dlon[dlon>180]-360) # convert to longitude ranging fromm -180:180
  
  CellAreaAtmos<-var.get.nc(ncin,"areacella",start=c(NA,NA),count=c(NA,NA))[maskLon,maskLat] # matrix rows(longitude) cols (latitude)
  CellA_bnds_lon<-var.get.nc(ncin,"lon_bnds",start=c(NA,NA),count=c(NA,NA))[,maskLon] # matrix rows(bounds) cols (longitude) - start 0 end 360
  CellA_bnds_lat<-var.get.nc(ncin,"lat_bnds",start=c(NA,NA),count=c(NA,NA))[,maskLat] # matrix rows(bounds) cols (latitude) - start minLat end maxLat
  
  # create dataframe with cell data in order for a raster and have bounds as well
  # lon, lat, value, lonmin,lonmax, latmin, latmax
  blat<-do.call(rbind,lapply(c(1:ncol(CellA_bnds_lat)),function(i,blat,n){matrix(rep(blat[,i],n),ncol=2,byrow=TRUE)},CellA_bnds_lat,nLon))
  blon<-matrix(CellA_bnds_lon,nrow=nLon*nLat,ncol=2,byrow=TRUE)
  blon[blon>180]<-(blon[blon>180]-360) # convert to longitude ranging fromm -180:180
  
  # dataframe from ACCESS atmosphere grid - vector of weights should correspond with rows in this grid
  dfA<-data.frame(lon=as.vector(dlon),lat=as.vector(dlat),area=as.vector(CellAreaAtmos),xmin=blon[,1],xmax=blon[,2],ymin=blat[,1],ymax=blat[,2])
  
  # note,
  #    dfA is not sorted by lat,lon so that lon orders from -178.125 to 180.  This should not affect outcome
  #    cells on 180 degrees have boundaries from 179.0625 to 180.9375 requiring a split of the cell at 180
  #    cells are created at -180 and appended as rows to account for these areas.  These are recombined in each MEASO area
  
  mask180<-dfA[,"lon"]==180
  dfA[mask180,"area"]<-dfA[mask180,"area"]/2  # split cell & will be rejoined later
  dfA_add<-dfA[mask180,]
  dfA_add[,c("lon","xmin")]<-c(-180,-180)
  dfA[mask180,"xmax"]<-180
  dfA<-rbind(dfA,dfA_add)
  
  # check if needed by converting to raster and accentuating shifts
  check<-FALSE
    if(check){
       tmp<-dfA[order(dfA[,"lat"],dfA[,"lon"]),]
       tmp[,"area"]<-tmp[,"area"]*1E-9
       tmp[tmp[,"lon"]<0,"area"]<-tmp[tmp[,"lon"]<0,"area"]*2
       r<-rast(tmp,type="xyz")
       crs(r)<-"epsg:4326"
       plot(r)
       } # end if
  
  # create polygons for cells 
 
    vCells<-vect(lapply(c(1:nrow(dfA)),function(i,dfA){
          p<-cbind(c(dfA[i,"xmin"],dfA[i,"xmin"],dfA[i,"xmax"],dfA[i,"xmax"],dfA[i,"xmin"])
             ,c(dfA[i,"ymin"],dfA[i,"ymax"],dfA[i,"ymax"],dfA[i,"ymin"],dfA[i,"ymin"]))
          return(vect(p,type="polygons",crs="epsg:4326"))
          },dfA))
    values(vCells)<-data.frame(area=dfA[,"area"],expanse=expanse(vCells),cell=seq(1,length(vCells),1))
  #    plot(vCells)
  
  MEASOshp<-vect(paste0(rootMEASOshp,fMEASOshp))
  
  # use if check needed
#  isRelated<-is.related(vCells,MEASOshp[1], "intersects")
#  plot(vCells[isRelated])
  
  MEASOcellWtsAtmos<-do.call(cbind,lapply(c(1:length(MEASOshp)),function(i,m,vC){
    isRelated<-is.related(vC,m[i], "intersects")
    ret<-rep(0,length(isRelated))
    iC<-intersect(m[i],vC)
    aC<-expanse(iC)/iC$expanse*iC$area
    ret[isRelated]<-aC[order(iC$cell)]
    ret[dfA[,"lon"]==180]<-ret[dfA[,"lon"]==180]+ret[dfA[,"lon"]==(-180)]
    return(ret[!(dfA[,"lon"]==(-180))])
  },MEASOshp,vCells))

  check<-FALSE
  if(check){
    tmp<-cbind(dfA[,"lon"],dfA[,"lat"],MEASOcellWtsAtmos[,15])
    r<-rast(tmp,type="xyz")
    crs(r)<-"epsg:4326"
    plot(r)
  } # end if

  return(MEASOcellWtsAtmos)  
} # end function - atmosphere grid


# 4. Generate Statistics ####
wtsO<-fnWtsOcean()
wtsA<-fnWtsAtmos()

MEASOsummary<-do.call(rbind,lapply(useVar,function(v,m,fV,wO,wA){
  root<-paste0(fV$root,fV[[v]]$var,"_y")
  ext<-paste0("_m",paste(sprintf("%02d",m),collapse="_"),fV$ext)
  vVar<-fV[[v]]$var
  vFiles<-fV[[v]][-1]
  nYrs<-length(vFiles)
  
  
  r1<-do.call(rbind,lapply(c(1:nYrs),function(i,vF,vVar,v,m,fV,wO,wA){
          y<-as.numeric(substr(names(vF)[i],2,5))
          dV<-do.call(cbind,readRDS(paste0(root,vF[[i]],ext)))

       if(v=="light"){ # median in space - average energy per day at the surface
            # use weights wA
             return(do.call(rbind,lapply(c(1:ncol(wA)),function(i,mCwts,dV,y){
               return(data.frame(var = v,measo = i,year = y,data = PCweighted(c(0.5),dV, mCwts[,i])))
            },wA,dV[,vVar],y)))
       } else if(v=="MLD") { # median in space - average mixed layer depth
             # use weights wO
             return(do.call(rbind,lapply(c(1:ncol(wO)),function(i,mCwts,dV,y){
               return(data.frame(var = v,measo = i,year = y,data = PCweighted(c(0.5),dV, mCwts[,i])))
            },wO,dV[,vVar],y)))
       } else if(v=="SeaIce") { # weighted mean - sea ice concentration (prop of area that is sea ice)
              # use weights wO
              return(do.call(rbind,lapply(c(1:ncol(wO)),function(i,mCwts,dV,y){
                return(data.frame(var = v,measo = i,year = y,data = MeanWeighted(dV, mCwts[,i])))
            },wO,dV[,vVar],y)))
       } else if(v=="temp") { # median in space - average temperature
              # use weights wO
              return(do.call(rbind,lapply(c(1:ncol(wO)),function(i,mCwts,dV,y){
                return(data.frame(var = v,measo = i,year = y,data = PCweighted(c(0.5),dV, mCwts[,i])))
            },wO,dV[,vVar],y)))
       } # end temp
  },vFiles,vVar,v,m,fV,wO,wA)) # end variable files
  return(r1)
  },useMonth,fVar,wtsO,wtsA))

saveRDS(MEASOsummary,"MEASOareaTimeSeriesSummary.rds")


# create nMDS matrix - standardise data sets
dnmds<-as.data.frame(do.call(cbind,lapply(useVar,function(v,d){
  d1<-d[d[,"var"]==v,"data"]
  print(d1)
  return((d1-mean(d1,na.rm=TRUE))/sd(d1,na.rm=TRUE))
  },MEASOsummary)))
names(dnmds)<-useVar

set.seed(123)
rnmds=metaMDS(dnmds,distance="euclidean")
rnmds
snmds<-as.data.frame(scores(rnmds))
snmds<-cbind(MEASOsummary[MEASOsummary[,"var"]==useVar[1],c("measo","year")],snmds)
snmds<-snmds[order(snmds[,"measo"],snmds[,"year"]),]
snmds[,"MEASOarea"]<-factor(MEASOnames[snmds[,"measo"]],levels = c("AOA", "AOS", "AON", "CIA", "CIS", "CIN", "EIA", "EIS", "EIN", "EPA", "EPS", "EPN", "WPA", "WPS", "WPN"))


plotPoints1900<-snmds[snmds[,"year"]==1900,]
plotPoints2020<-snmds[snmds[,"year"]==2020,]

# MEASO colours ####
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

pointCols<-as.vector(c(MEASOcols[2],MEASOcols[5],MEASOcols[8],MEASOcols[11],MEASOcols[14]))
names(pointCols)<-NULL #c("AO","CI","EI","EP","WP")

SectorColours<-pointCols
MEASOcolours<-do.call(c,lapply(SectorColours,function(sc){rep(sc,3)}))
MEASOpoints<-rep(c(1,0,5),5)


p <- ggplot(snmds, aes(x = NMDS1, y = NMDS2)) 
p <- p+  geom_path(aes(colour = MEASOarea))
p <- p+ scale_color_manual(values=MEASOcolours)
p <- p+ geom_point(stat="identity",aes(x=NMDS1,y=NMDS2,colour=MEASOarea,fill=MEASOarea),size=1)
p <- p+ geom_point(data=plotPoints1900
                   ,stat="identity",aes(x=NMDS1,y=NMDS2,colour=MEASOarea,fill=MEASOarea,shape=MEASOarea),size=5)
p <- p+ geom_point(data=plotPoints2020
                   ,stat="identity",aes(x=NMDS1,y=NMDS2,colour=MEASOarea,fill=MEASOarea,shape=MEASOarea),size=3)
p <- p+ scale_shape_manual(values=MEASOpoints) #,labels=MEASOnames)
#p <- p+ geom_text(x=1, y=2, label="Stress < 0.001")
p <- p+ theme(axis.text.y = element_text(colour = "black", size = 12)
              , axis.text.x = element_text(colour = "black", size = 12)
              ,legend.text = element_text(size = 12, colour ="black") 
              ,legend.position = "right"
              ,legend.title = element_text( size = 14, colour = "black")
#              ,axis.title.y = element_text(face = "bold", size = 14) 
#              ,axis.title.x = element_text(face = "bold", size = 14, colour = "black") 
              ,panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2)
              ,legend.key=element_blank()
              )# end theme
p <- p+ labs(x = "", y = "")
p

