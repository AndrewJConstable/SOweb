# ACCESS ESM1.5 SSP5-8.5
# supported by Tilo Ziehn

# Variable definitions:
# https://clipc-services.ceda.ac.uk/dreq/mipVars.html

#October 2021: Tilo provide surface temperature and sea ice monthly means for COP 26 animation
#              also provided global warming level (see emails from 26 October 2021)

# https://dapds00.nci.org.au/thredds/catalog/fs38/publications/CMIP6/C4MIP/CSIRO/ACCESS-ESM1-5/esm-ssp585/r1i1p1f1/Omon/catalog.html

library(ncdf4)
library(terra)

ncin <- nc_open("/Volumes/WizardAtWork/_d/ACCESS\ ESM/SO_1900-2100/SIC/MonthlyMeans_SIC_1850-2014.nc")
print(ncin)

