# routine to process files from
# Fons S, Kurtz N, Bagnardi M. A decade-plus of Antarctic sea ice thickness and volume 
# estimates from CryoSat-2 using a physical model and waveform fitting. 
# The Cryosphere. 2023;17(6):2487-508. (https://tc.copernicus.org/articles/17/2487/2023/)
# obtained from data source: https://zenodo.org/records/7327711
# 
library(RNetCDF)
library(terra)


fPath<-"/Users/acon/Desktop/_w/_d/Ocean/SeaIce/Fons et al 2023/"
fName<-"CS2WFA_25km_201007.nc"

tmp1<-rast(paste0(fPath,fName))
tmp2<-subset(tmp1,"sea_ice_thickness")

# files have no CRS, resolution or extent (imported as netCDF or raster)
# notes say netCDF based on NSIDC polar stereographic 25km resolution
# retrieved projection and extent from NSIDC website (https://nsidc.org/data/user-resources/help-center/guide-nsidcs-polar-stereographic-projection)
crs(tmp2)<-"EPSG:3976"
ext(tmp2)<-c(-3950,3950,-3950,4350) # obtained from NSIDC website on their polar stereographic grid
values(tmp2)[is.nan(values(tmp2))]<-NA

# code thanks to Mike Sumner 
r <- flip(rast(paste0(fPath,fName)), "vertical")
set.ext(r, ext(-3950000, 3950000, -3950000, 4350000))
set.crs(r, "EPSG:3031")
targetgrid <- rast(ext(-180, 180, -90, -50), res = c(0.25, 0.25), crs = "EPSG:4326")
plot(project(r, targetgrid, by_util = TRUE))

plot(project(subset(r,"sea_ice_thickness"), targetgrid, by_util = TRUE))

