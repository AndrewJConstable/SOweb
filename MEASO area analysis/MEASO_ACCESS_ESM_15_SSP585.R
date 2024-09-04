# based on MEASO BRAN2022 Temp means - grp[yr] grp[mth] grp[depth].R

# function to generate raster of mean temperatures across a group of years, months and depths
# in following sequence
#. 1. weighted mean temperature between min and max depth (can use mixed layer depth as either minimum or maximum)
#  2. mean temperature across months as designated in input vector
#. 3. mean temperature across years as designated in input vector

# Important notes for processing BRAN data
# 1. Using Terra package to import netCDF files.
# 2. Rasters are complete x-y matrices, including when subsetting
# 3. When converting rasters to dataframes using as.data.frame, missing values are dropped if the defaults are used.
#          Important to not remove cells with NA by setting na.rm=F i.e. as.data.frame(r, na.rm=F)
#          and manage NAs in code.

# Important notes for processing ACCESS ESM1.5 SSP585 replicate 1

#########################################################################################################
# 1. Libraries
##################################

library(lubridate)
library(RNetCDF)
#library(ncdf4) # for reading the netCDF files which are on an irregular grid
library(interp)
library(terra)



# calculating pH, pCO2 from Temp, alkalinity, DIC
# use seacarb library  (utube tutorial : https://www.youtube.com/watch?v=35_Utg8hiAM&t=86s)

# NCI data: CMIP6 ACCESS ESM 1.5
# historical runs 1850-2014: https://dapds00.nci.org.au/thredds/catalog/fs38/publications/CMIP6/CMIP/CSIRO/ACCESS-ESM1-5/esm-hist/catalog.html
# SSP585 2015-2100: https://dapds00.nci.org.au/thredds/catalog/fs38/publications/CMIP6/ScenarioMIP/CSIRO/ACCESS-ESM1-5/catalog.html


# light - shortwave flux at surface (W.m-2)

esmLight<-list(files=c( "https://dapds00.nci.org.au/thredds/dodsC/fs38/publications/CMIP6/CMIP/CSIRO/ACCESS-ESM1-5/esm-hist/r1i1p1f1/Amon/rsds/gn/latest/rsds_Amon_ACCESS-ESM1-5_esm-hist_r1i1p1f1_gn_185001-201412.nc"
          ,"https://dapds00.nci.org.au/thredds/dodsC/fs38/publications/CMIP6/C4MIP/CSIRO/ACCESS-ESM1-5/esm-ssp585/r1i1p1f1/Amon/rsds/gn/latest/rsds_Amon_ACCESS-ESM1-5_esm-ssp585_r1i1p1f1_gn_201501-210012.nc")
         ,var="rsds")

esmIce <-list(files=c("https://dapds00.nci.org.au/thredds/dodsC/fs38/publications/CMIP6/CMIP/CSIRO/ACCESS-ESM1-5/esm-hist/r1i1p1f1/SImon/siconc/gn/latest/siconc_SImon_ACCESS-ESM1-5_esm-hist_r1i1p1f1_gn_185001-201412.nc"
           ,"https://dapds00.nci.org.au/thredds/dodsC/fs38/publications/CMIP6/C4MIP/CSIRO/ACCESS-ESM1-5/esm-ssp585/r1i1p1f1/SImon/siconc/gn/latest/siconc_SImon_ACCESS-ESM1-5_esm-ssp585_r1i1p1f1_gn_201501-210012.nc")
         ,var = "siconc")

esmTemp <- list(files=c( "https://dapds00.nci.org.au/thredds/dodsC/fs38/publications/CMIP6/CMIP/CSIRO/ACCESS-ESM1-5/esm-hist/r10i1p1f1/Omon/thetao/gn/latest/thetao_Omon_ACCESS-ESM1-5_esm-hist_r10i1p1f1_gn_185001-185912.nc"
           ,"https://dapds00.nci.org.au/thredds/dodsC/fs38/publications/CMIP6/CMIP/CSIRO/ACCESS-ESM1-5/esm-hist/r10i1p1f1/Omon/thetao/gn/latest/thetao_Omon_ACCESS-ESM1-5_esm-hist_r10i1p1f1_gn_186001-186912.nc"
           ,"https://dapds00.nci.org.au/thredds/dodsC/fs38/publications/CMIP6/CMIP/CSIRO/ACCESS-ESM1-5/esm-hist/r10i1p1f1/Omon/thetao/gn/latest/thetao_Omon_ACCESS-ESM1-5_esm-hist_r10i1p1f1_gn_187001-187912.nc"
           ,"https://dapds00.nci.org.au/thredds/dodsC/fs38/publications/CMIP6/CMIP/CSIRO/ACCESS-ESM1-5/esm-hist/r10i1p1f1/Omon/thetao/gn/latest/thetao_Omon_ACCESS-ESM1-5_esm-hist_r10i1p1f1_gn_188001-188912.nc"
           ,"https://dapds00.nci.org.au/thredds/dodsC/fs38/publications/CMIP6/CMIP/CSIRO/ACCESS-ESM1-5/esm-hist/r10i1p1f1/Omon/thetao/gn/latest/thetao_Omon_ACCESS-ESM1-5_esm-hist_r10i1p1f1_gn_189001-189912.nc"
           ,"https://dapds00.nci.org.au/thredds/dodsC/fs38/publications/CMIP6/CMIP/CSIRO/ACCESS-ESM1-5/esm-hist/r10i1p1f1/Omon/thetao/gn/latest/thetao_Omon_ACCESS-ESM1-5_esm-hist_r10i1p1f1_gn_190001-190912.nc"
           ,"https://dapds00.nci.org.au/thredds/dodsC/fs38/publications/CMIP6/CMIP/CSIRO/ACCESS-ESM1-5/esm-hist/r10i1p1f1/Omon/thetao/gn/latest/thetao_Omon_ACCESS-ESM1-5_esm-hist_r10i1p1f1_gn_191001-191912.nc"
           ,"https://dapds00.nci.org.au/thredds/dodsC/fs38/publications/CMIP6/CMIP/CSIRO/ACCESS-ESM1-5/esm-hist/r10i1p1f1/Omon/thetao/gn/latest/thetao_Omon_ACCESS-ESM1-5_esm-hist_r10i1p1f1_gn_192001-192912.nc"
           ,"https://dapds00.nci.org.au/thredds/dodsC/fs38/publications/CMIP6/CMIP/CSIRO/ACCESS-ESM1-5/esm-hist/r10i1p1f1/Omon/thetao/gn/latest/thetao_Omon_ACCESS-ESM1-5_esm-hist_r10i1p1f1_gn_193001-193912.nc"
           ,"https://dapds00.nci.org.au/thredds/dodsC/fs38/publications/CMIP6/CMIP/CSIRO/ACCESS-ESM1-5/esm-hist/r10i1p1f1/Omon/thetao/gn/latest/thetao_Omon_ACCESS-ESM1-5_esm-hist_r10i1p1f1_gn_194001-194912.nc"
           ,"https://dapds00.nci.org.au/thredds/dodsC/fs38/publications/CMIP6/CMIP/CSIRO/ACCESS-ESM1-5/esm-hist/r10i1p1f1/Omon/thetao/gn/latest/thetao_Omon_ACCESS-ESM1-5_esm-hist_r10i1p1f1_gn_195001-195912.nc"
           ,"https://dapds00.nci.org.au/thredds/dodsC/fs38/publications/CMIP6/CMIP/CSIRO/ACCESS-ESM1-5/esm-hist/r10i1p1f1/Omon/thetao/gn/latest/thetao_Omon_ACCESS-ESM1-5_esm-hist_r10i1p1f1_gn_196001-196912.nc"
           ,"https://dapds00.nci.org.au/thredds/dodsC/fs38/publications/CMIP6/CMIP/CSIRO/ACCESS-ESM1-5/esm-hist/r10i1p1f1/Omon/thetao/gn/latest/thetao_Omon_ACCESS-ESM1-5_esm-hist_r10i1p1f1_gn_197001-197912.nc"
           ,"https://dapds00.nci.org.au/thredds/dodsC/fs38/publications/CMIP6/CMIP/CSIRO/ACCESS-ESM1-5/esm-hist/r10i1p1f1/Omon/thetao/gn/latest/thetao_Omon_ACCESS-ESM1-5_esm-hist_r10i1p1f1_gn_198001-198912.nc"
           ,"https://dapds00.nci.org.au/thredds/dodsC/fs38/publications/CMIP6/CMIP/CSIRO/ACCESS-ESM1-5/esm-hist/r10i1p1f1/Omon/thetao/gn/latest/thetao_Omon_ACCESS-ESM1-5_esm-hist_r10i1p1f1_gn_199001-199912.nc"
           ,"https://dapds00.nci.org.au/thredds/dodsC/fs38/publications/CMIP6/CMIP/CSIRO/ACCESS-ESM1-5/esm-hist/r10i1p1f1/Omon/thetao/gn/latest/thetao_Omon_ACCESS-ESM1-5_esm-hist_r10i1p1f1_gn_200001-200912.nc"
           ,"https://dapds00.nci.org.au/thredds/dodsC/fs38/publications/CMIP6/CMIP/CSIRO/ACCESS-ESM1-5/esm-hist/r10i1p1f1/Omon/thetao/gn/latest/thetao_Omon_ACCESS-ESM1-5_esm-hist_r10i1p1f1_gn_201001-201412.nc"
           ,"https://dapds00.nci.org.au/thredds/dodsC/fs38/publications/CMIP6/C4MIP/CSIRO/ACCESS-ESM1-5/esm-ssp585/r1i1p1f1/Omon/thetao/gn/latest/thetao_Omon_ACCESS-ESM1-5_esm-ssp585_r1i1p1f1_gn_201501-202412.nc"
           ,"https://dapds00.nci.org.au/thredds/dodsC/fs38/publications/CMIP6/C4MIP/CSIRO/ACCESS-ESM1-5/esm-ssp585/r1i1p1f1/Omon/thetao/gn/latest/thetao_Omon_ACCESS-ESM1-5_esm-ssp585_r1i1p1f1_gn_202501-203412.nc"
           ,"https://dapds00.nci.org.au/thredds/dodsC/fs38/publications/CMIP6/C4MIP/CSIRO/ACCESS-ESM1-5/esm-ssp585/r1i1p1f1/Omon/thetao/gn/latest/thetao_Omon_ACCESS-ESM1-5_esm-ssp585_r1i1p1f1_gn_203501-204412.nc"
           ,"https://dapds00.nci.org.au/thredds/dodsC/fs38/publications/CMIP6/C4MIP/CSIRO/ACCESS-ESM1-5/esm-ssp585/r1i1p1f1/Omon/thetao/gn/latest/thetao_Omon_ACCESS-ESM1-5_esm-ssp585_r1i1p1f1_gn_204501-205412.nc"
           ,"https://dapds00.nci.org.au/thredds/dodsC/fs38/publications/CMIP6/C4MIP/CSIRO/ACCESS-ESM1-5/esm-ssp585/r1i1p1f1/Omon/thetao/gn/latest/thetao_Omon_ACCESS-ESM1-5_esm-ssp585_r1i1p1f1_gn_205501-206412.nc"
           ,"https://dapds00.nci.org.au/thredds/dodsC/fs38/publications/CMIP6/C4MIP/CSIRO/ACCESS-ESM1-5/esm-ssp585/r1i1p1f1/Omon/thetao/gn/latest/thetao_Omon_ACCESS-ESM1-5_esm-ssp585_r1i1p1f1_gn_206501-207412.nc"
           ,"https://dapds00.nci.org.au/thredds/dodsC/fs38/publications/CMIP6/C4MIP/CSIRO/ACCESS-ESM1-5/esm-ssp585/r1i1p1f1/Omon/thetao/gn/latest/thetao_Omon_ACCESS-ESM1-5_esm-ssp585_r1i1p1f1_gn_207501-208412.nc"
           ,"https://dapds00.nci.org.au/thredds/dodsC/fs38/publications/CMIP6/C4MIP/CSIRO/ACCESS-ESM1-5/esm-ssp585/r1i1p1f1/Omon/thetao/gn/latest/thetao_Omon_ACCESS-ESM1-5_esm-ssp585_r1i1p1f1_gn_208501-209412.nc"
           ,"https://dapds00.nci.org.au/thredds/dodsC/fs38/publications/CMIP6/C4MIP/CSIRO/ACCESS-ESM1-5/esm-ssp585/r1i1p1f1/Omon/thetao/gn/latest/thetao_Omon_ACCESS-ESM1-5_esm-ssp585_r1i1p1f1_gn_209501-210012.nc")
          ,var="thetao")

esmMLD <- list(files = c(  "https://dapds00.nci.org.au/thredds/dodsC/fs38/publications/CMIP6/CMIP/CSIRO/ACCESS-ESM1-5/esm-hist/r10i1p1f1/Omon/mlotst/gn/latest/mlotst_Omon_ACCESS-ESM1-5_esm-hist_r10i1p1f1_gn_185001-201412.nc"
           ,"https://dapds00.nci.org.au/thredds/dodsC/fs38/publications/CMIP6/C4MIP/CSIRO/ACCESS-ESM1-5/esm-ssp585/r1i1p1f1/Omon/mlotst/gn/latest/mlotst_Omon_ACCESS-ESM1-5_esm-ssp585_r1i1p1f1_gn_201501-210012.nc")
          ,var = "mlotst")


#########################################################################################################
# 3. Functions
##################################
domainBounds<-list(limLon=c(0,360),limLat=c(-80,-30),resolution=0.1) # also for regular grid to use

inputYrs<-list(early=c(1981:2000),late=c(2080:2099))
inputMths<-list(period1=c(10),period2=c(11),period3=c(12),period4=c(10,11),period5=c(11,12)) # order counts - if january is after december then january is drawn from the next year

inputVar<-list(esmLight,esmIce,esmTemp,esmMLD)

for(i in c(1:length(inputYrs))){
  Yrs<-inputYrs[[i]]
  for(j in c(1:length(inputMths))){
    Mths<-inputMths[[j]]
    for(k in c(1:length(inputVar))){
      fESM<-inputVar[[k]]$files
      vESM<-inputVar[[k]]$var
    print(paste(Yrs[1],Mths[1],vESM))
    summariseESMvar(vESM,fESM,Yrs,Mths,domainBounds)
}}}  


tmp<-summariseESMvar(vESM,fESM,Yrs,Mths,domainBounds)

############# start main function #######################
summariseESMvar<-function(vESM,fESM,Yrs,Mths,domBds){

# 1. loop through files and generate times in files
  print("1")
  
fileDates<-do.call(rbind,lapply(c(1:length(fESM)),function(i,f){
      ncin<-open.nc(f[i])
      time<-var.get.nc(ncin,"time")
      t<-as_date(time,"1850-01-01")
      Y<-year(t)
      M<-month(t)
      fi<-rep(i,length(t))
      df<-data.frame(Date=t,Year=Y,Month=M,FileNum=fi,Time=c(1:length(t)))
      return(df)
  },fESM)) # end fileDates

# 2. generate year list of file-time combinations to be averaged
print("2")
detYrMth<-lapply(Yrs,function(y,m,f,fileDates){
    # create matrix of year,month to allow for months in following year e.g. c(11,12,1,4)
  nMths<-length(m)
  YrMth<-matrix(c(rep(y,nMths),m),ncol=2,dimnames=list(NULL,c("tgtYr","month")))
  if(nMths>1){
      doNextYr<-m[c(2:nMths)]<m[c(1:(nMths-1))]
      if(sum(doNextYr)>0){
       doNextYr<-which(doNextYr)[1]+1
       YrMth[c(doNextYr:nMths),"tgtYr"]<-y+1
       }
    } # end nMths>1
  return(matrix(unlist(apply(YrMth,1,function(ym,fd){fd[which(fd[,"Year"]==ym[1] & fd[,"Month"]==ym[2]),c("FileNum","Time")]},fileDates)),ncol=2,byrow=TRUE,dimnames=list(NULL,c("FileNum","Time"))))
    },Mths,fESM,fileDates)

# 3. from first file and first time, generate domain mask and longitudes and latitudes for final dataframe
print("3")

domain<-{
     ncin<-open.nc(fESM[1])
     if(vESM=="rsds") { # regular grid
       lat<-var.get.nc(ncin,"lat",start=c(NA),count=c(NA))
       lon<-var.get.nc(ncin,"lon",start=c(NA),count=c(NA))
       nlat<-length(lat)
       nlon<-length(lon)
       lat<-matrix(lat,nrow=nlon,ncol=nlat,byrow=TRUE)
       lon<-matrix(lon,nrow=nlon,ncol=nlat)
     } else { # irregular grid
         lat<-var.get.nc(ncin,"latitude",start=c(NA,NA),count=c(NA,NA))
         lon<-var.get.nc(ncin,"longitude",start=c(NA,NA),count=c(NA,NA))
         }
     close.nc(ncin)
     maskDomain<- lat>=domBds$limLat[1] & lat<=domBds$limLat[2] & lon>=domBds$limLon[1] & lon<=domBds$limLon[2]
     list(lon=lon[maskDomain],lat=lat[maskDomain],mask=maskDomain)
      }
# 4. loop through years and average the months
#   note that the domain mask is sufficient for generating columns for each month - then take mean of rows.
print("4")
resYrs <- sapply(detYrMth,function(dym,vESM,fESM,limLon,limLat,resolution,maskDomain){

     # dym  matrix - rows(month) cols(Filenum,Time)

     res<-apply(dym,1,function(fT,fESM,mask,vESM){
            ncin<-open.nc(fESM[fT[1]])
            # check if depth is part of data array
            tmp1<-read.nc(ncin,start=rep(1,10),count=rep(1,10))
            if(sum(names(tmp1) %in% "lev")>0){
            dVar<-var.get.nc(ncin,vESM,start=c(NA,NA,1,fT[2]),count=c(NA,NA,1,1)) # (i,j,lev,time)
            } else {
            dVar<-var.get.nc(ncin,vESM,start=c(NA,NA,fT[2]),count=c(NA,NA,1)) # (i,j,time)  
            }
            close.nc(ncin)
            return(dVar[maskDomain])
                },fESM,maskDomain,vESM)
     return(apply(res,1,mean))
    },vESM,fESM,domBds$limLon,domBds$limLat,resolution,domain$mask)

resYM<-list(domain$lon,domain$lat,apply(resYrs,1,mean))
names(resYM)<-c("lon","lat",vESM)
saveRDS(resYM,file=paste0(vESM,"_y",Yrs[1],"-",Yrs[length(Yrs)],"_m",paste(sprintf("%02d",Mths),sep="",collapse="_"),".Rdata"))
return(resYM)
} # end function






dVar<-tmp$rsds
dLon<-tmp$lon
dLat<-tmp$lat

# check
# prepare grid
gridLat<-seq((domainBounds$limLat[1]+domainBounds$resolution/2),(domainBounds$limLat[2]-domainBounds$resolution/2),domainBounds$resolution)
gridLon<-seq((domainBounds$limLon[1]+domainBounds$resolution/2),(domainBounds$limLon[2]-domainBounds$resolution/2),domainBounds$resolution)

# create gridded map of missing cells
maskMissing<-is.na(dVar)
dVarMiss<-dVar *0+1
dVarMiss[maskMissing]<-0
gridMiss<-interp(x=dLon,y=dLat,z=dVarMiss,xo=gridLon,yo=gridLat,method="linear",na.rm=TRUE)

gdfMiss<-data.frame(lon=as.vector(matrix(gridMiss$x,nrow=length(gridLon),ncol=length(gridLat)))
                ,lat=as.vector(matrix(gridMiss$y,nrow=length(gridLon),ncol=length(gridLat),byrow=TRUE))
                ,miss=as.vector(gridMiss$z))
rastMiss<-rast(gdfMiss)
rastMiss[rastMiss<0.7]<-NA
rastMiss[!is.na(rastMiss)]<-1
plot(rastMiss)

# grid data
mask<-!maskMissing

gridData<-interp(x=dLon[mask],y=dLat[mask],z=dVar[mask],xo=gridLon,yo=gridLat,na.rm=TRUE)

gdf<-data.frame(lon=as.vector(matrix(gridData$x,nrow=length(gridLon),ncol=length(gridLat)))
                ,lat=as.vector(matrix(gridData$y,nrow=length(gridLon),ncol=length(gridLat),byrow=TRUE))
                ,Var=as.vector(gridData$z))
rastRes<-rast(gdf)
plot(rastRes)
plot(rastRes*rastMiss)


