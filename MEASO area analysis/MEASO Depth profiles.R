# Process and save depth profiles for MEASO areas
library(ggplot2)
library(ggforce)
library(terra)

# MEASO areas
rootMEASOshp<-"/Users/acon/Desktop/_w/_d/Shapefiles/MEASO polygons/"
fMEASOshp<-"MEASO_polygons.shp"

MEASOshp<-vect(paste0(rootMEASOshp,fMEASOshp))

# GEBCO data 2023

gebco23<-rast("/Users/acon/Desktop/_w/_d/GEBCO_10_Apr_2024_41650c45d0e9/GEBCO23 Southern Ocean.tif")
names(gebco23)<-"height"
gebco23cellarea<-cellSize(gebco23,unit="km")
gebco23<-c(gebco23,gebco23cellarea)
#plot(gebco23)

# GEBCO23 data - read raw tif files (could only download as tiles because single download did not work), merge and write
# gebco23_1<-rast("/Users/acon/Desktop/_w/_d/GEBCO_10_Apr_2024_41650c45d0e9/gebco_2023_n-30.0_s-90.0_w-180.0_e-60.0.tif")
# gebco23_2<-rast("/Users/acon/Desktop/_w/_d/GEBCO_10_Apr_2024_41650c45d0e9/gebco_2023_n-30.0_s-90.0_w-65.0_e5.0.tif")
# gebco23_3<-rast("/Users/acon/Desktop/_w/_d/GEBCO_10_Apr_2024_41650c45d0e9/gebco_2023_n-30.0_s-90.0_w0.0_e65.0.tif")
# gebco23_4<-rast("/Users/acon/Desktop/_w/_d/GEBCO_10_Apr_2024_41650c45d0e9/gebco_2023_n-30.0_s-90.0_w60.0_e180.0.tif")
# gebco23<-merge(gebco23_1,gebco23_2)
# gebco23<-merge(gebco23,gebco23_3)
# gebco23<-merge(gebco23,gebco23_4)
# writeRaster(gebco23,"GEBCO23.tif") # moved to data directory

MEASOshp<-project(MEASOshp,gebco23)
plot(gebco23$height)
plot(MEASOshp,add=TRUE)

depthInt<-c(-1000,-2000,-3000,-4000)

dDepth<-lapply(seq(1,length(MEASOshp),1),function(m,ma,g,di){
  Ageom<-geom(ma[m])
  Nparts<-unique(Ageom[,"part"]) # separate of polygon parts (around -180) is important because of limits in vector processing
  mDepth<-do.call(rbind,lapply(c(1:length(Nparts)),function(p,g,a){
     a1<-a[a[,"part"]==p,]
     a1[,"part"]<-1
     a1<-vect(a1,type="polygons",crs=crs(g))
     g1<-crop(g,a1)
     return(extract(g1,a1,cells=TRUE,weights=TRUE))     
     },g,Ageom))
  
  CellWt <- mDepth[,"area"]*mDepth[,"weight"]
  HeightMax<-max(mDepth[,"height"],na.rm=TRUE)
  A      <- sum(CellWt,na.rm=TRUE)
  binA   <-sapply(c(0,di),function(di,cw,d){
    useD   <- d >= di
    return(sum(cw[useD]))
  },CellWt,mDepth[,"height"])
  return(list(TotalArea = A, MaxHeight = HeightMax, depth=c(0,di), prop = binA/A))
},MEASOshp,gebco23,depthInt)

saveRDS(dDepth,file="DepthProfile.rds")
