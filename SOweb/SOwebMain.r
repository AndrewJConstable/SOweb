# 0. Setup ####

#    0.1 Logic of variable construction ####

#         Constants = input data held constant
#         Temporal vector = input data used in integration.  
#                           Element position relates to the value applied during a time step
#         Matrix of State Variables = rows(State Variables), cols(End of time step)
#                                     ncols = nTimeSteps where the last column is the end
#                                             of the integration period
#         Initial State Vector = length matches the number of rows in the matrix above
#         

#    0.2 Libraries ####

library(ggplot2)
library(gplots)
library(spatstat)
library(stats4)
library(bbmle)

TargetArea <- "AOA" # CIS

SourceFiles<-list( Fns_General       = "SOwebFns_General.r"
                  ,Fns_Detritus      = "SOwebFns_Detritus.r"
                  ,Fns_Phytoplankton = "SOwebFns_Phytoplankton.r"
                  ,Fns_Nutrients     = "SOwebFns_Nutrients.r"
                  ,Data_cE           = "SOwebData_cE.r"
                  ,Data_pools        = "SOwebData_pools.r"
   )# end source files

# 1. Minimisation inputs

# time0 relates to 1st element of a temporal vector.  The matrix of state variables 
nYrs<-1
nStepsPerYr<-365 
Date1<-as.Date("2015-10-01")  # start date chosen to be at maximum sea ice when detritus likely to be at a minimum

# 2.1 Managing time ####
DaysInMonth<-c(31,28,31,30,31,30,31,31,30,31,30,31)
DayMidMonth<-c(15,14,15,15,15,15,15,15,15,15,15,15)


#    1.4. Input Environmental Data ####

fPathIn<-"./DataIn/"
fPathOut<-"./DataOut/"
load(paste0(fPathIn,"Fe_means.Rdata")) # df Iron concentration (micromol.m-3) rows(MEASO area) cols(n, name, value)
load(paste0(fPathIn,"Si_means.Rdata")) # df Silicate concentration (mmol.m-3) rows(MEASO area) cols(n, name, value)
load(paste0(fPathIn,"Lat_means.Rdata")) # df Latitude mean for MEASO areas rows(MEASO area) cols(n, name, value)
load(paste0(fPathIn,"Insol_means.Rdata"))  # matrix Insolation (W.m-2.d-1) rows(months) cols(MEASO area)
load(paste0(fPathIn,"MLD_means.Rdata"))  # matrix Mixed Layer Depth (m) rows(months) cols(MEASO area)
load(paste0(fPathIn,"CICE_means.Rdata"))  # matrix Sea ice percent cover (%) rows(months) cols(MEASO area) 
load(paste0(fPathIn,"Temp_MLD_means.Rdata"))  # matrix Temperature mixed layer (oC) rows(months) cols(MEASO area) 
load(paste0(fPathIn,"Temp_Deep_means.Rdata"))  # matrix Temperature MLD to 1000m (oC) rows(months) cols(MEASO area) 


# 3. Subset Environmental Data ####

nFe_conc_0      <-Fe_means[Fe_means[,"MEASO_area"]==TargetArea,3]   # micromole.m-3
nSi_conc_0      <-Si_means[Si_means[,"MEASO_area"]==TargetArea,3]   # millimole.m-3
dCaSI_conc_0    <- 0                                                # g.m-3
dCaM_conc_0     <- 0                                                # g.m-3
dCaD_conc_0     <- 0                                                # g.m-3
dCaS_conc_0     <- 0                                                # g.m-3
eLat       <-Lat_means[Lat_means[,"MEASO_area"]==TargetArea,3] # degrees (negative = south)
eInsol     <-Insol_means[,TargetArea]                          # Watts.m-2.day-1
eMLD       <-MLD_means[,TargetArea]                            # m
eCICE      <-CICE_means[,TargetArea]                           # percent
eTemp_MLD  <-Temp_MLD_means[,TargetArea]                       # oC
eTemp_Deep <-Temp_Deep_means[,TargetArea]                      # oC

# 2. Input functions and data ####

source(SourceFiles$Fns_General)
source(SourceFiles$Fns_Detritus)
source(SourceFiles$Fns_Phytoplankton)
source(SourceFiles$Fns_Nutrients)
source(SourceFiles$Data_cE)

cE<-c(cE,list(Lat = eLat
              ,Conc_t0 = list(
                           nFeM   = nFe_conc_0   # initial iron concentration umole.m-3 in mixed layer 
                          ,nFeD   = nFe_conc_0   # initial iron concentration umole.m-3 in deep
                          ,nFeSI  = nFe_conc_0   # initial iron concentration umole.m-3 in sea ice
                          ,nSiM   = nSi_conc_0   # initial silicic acid concentration mmole.m-3 in mixed layer
                          ,nSiD   = nSi_conc_0   # initial silicic acid concentration mmole.m-3 in deep
                          ,nSiSI  = nSi_conc_0   # initial silicic acid concentration mmole.m-3 in sea ice
                          ,dCaSI  = dCaSI_conc_0 # initial carbon concentration g.m-3 in sea ice
                          ,dCaM   = dCaM_conc_0  # initial carbon concentration g.m-3 in mixed layer
                          ,dCaD   = dCaD_conc_0  # initial carbon concentration g.m-3 in deep
                          ,dCaS   = dCaS_conc_0  # initial carbon concentration g.m-3 in sediment
                          ,dFeSI  = dCaSI_conc_0*cE$C_MassToMole*0.005 # initial dFe concentration umole.m-3 in sea ice (0.005 as qFe from Hauck et al 2013)
                          ,dFeM   = dCaM_conc_0*cE$C_MassToMole*0.005  # initial dFe concentration umole.m-3 in mixed layer
                          ,dFeD   = dCaD_conc_0*cE$C_MassToMole*0.005  # initial dFe concentration umole.m-3 in deep
                          ,dFeS   = dCaS_conc_0*cE$C_MassToMole*0.005  # initial dFe concentration umole.m-3 in sediment
                          ,dSiSI  = dCaSI_conc_0*cE$C_MassToMole*0.04 # initial dSi concentration mmole.m-3 in sea ice (0.04 as qSimin from Hauck et al 2013)
                          ,dSiM   = dCaM_conc_0*cE$C_MassToMole*0.04  # initial dSi concentration mmole.m-3 in mixed layer
                          ,dSiD   = dCaD_conc_0*cE$C_MassToMole*0.04  # initial dSi concentration mmole.m-3 in deep
                          ,dSiS   = dCaS_conc_0*cE$C_MassToMole*0.04  # initial dSi concentration mmole.m-3 in sediment
                          ) # end Conc_t0
) # end new list
      ) # end add to cE


source("SOwebData_pools.r") # here for using parameters in cE


# 2. Pre-minimisation simplification and standardisation ####

#    2.1 Managing time ####
day<-c(1:nStepsPerYr)
nTimeSteps<-nYrs*nStepsPerYr

#    2.2 Environmental variables ####
#
#        2.2.1 constants

# update cE

poolDepth<-sapply(names(a) # read depth layers of pools
                  ,function(s,a){if(is.null(a[[s]]$Attr$Layer)) NA else a[[s]]$Attr$Layer},a)

cE<-c(cE,list(Conc_t0_layer = poolDepth))

#        2.2.2 changing order of environmental variables to start with the first month ####
Month1<-format(Date1,"%m")

fnReorderMonthVector<-function(m,V){
       e<-which(names(eTemp_MLD)%in%paste0("M",m))
       if(e==1) return(V) else return(V[c(e:12,1:(e-1))]) 
       } # end fn
orderMonths<-fnReorderMonthVector(Month1,seq(1,12,1))
eInsol     <-fnReorderMonthVector(Month1,eInsol)
eMLD       <-fnReorderMonthVector(Month1,eMLD)
eCICE      <-fnReorderMonthVector(Month1,eCICE)
eTemp_MLD  <-fnReorderMonthVector(Month1,eTemp_MLD)
eTemp_Deep <-fnReorderMonthVector(Month1,eTemp_Deep)

#        2.2.3 interpolation of environmental vectors to time steps ####
# interpolation to all time steps based on the monthly values being mid month
# do one year and then replicate for the number of years needed

# vector for days corresponding to months plus beginning and end of year for interpolation
DayMM<-c(0,DayMidMonth[orderMonths][1],(cumsum(DaysInMonth[orderMonths])[1:11]+DayMidMonth[orderMonths][2:12]),366)

# vector of day in year at each time step
library(lubridate)
Day1<-yday(Date1)
detach(package:lubridate,unload=TRUE)
Jday<-seq(Day1,365,1)
if(Day1>1) Jday<-c(Jday,seq(1,(Day1-1),1))

tE<-list(day      = rep(Jday,nYrs)  # the time steps being used
        ,Insol     = rep(fnInterpolateVectorToTimeSteps(eInsol,DayMM,day),nYrs)
        ,MLD       = rep(fnInterpolateVectorToTimeSteps(eMLD,DayMM,day),nYrs)
        ,CICE      = rep(fnInterpolateVectorToTimeSteps(eCICE,DayMM,day),nYrs)/100 # as proportion rather than percent
        ,OpenWater = (1-rep(fnInterpolateVectorToTimeSteps(eCICE,DayMM,day),nYrs)/100) # as proportion rather than percent
        ,MLD       = rep(fnInterpolateVectorToTimeSteps(eTemp_MLD,DayMM,day),nYrs)
        ,Temp_MLD  = rep(fnInterpolateVectorToTimeSteps(eTemp_MLD,DayMM,day),nYrs)
        ,Temp_Deep = rep(fnInterpolateVectorToTimeSteps(eTemp_Deep,DayMM,day),nYrs)
         ) # end tE

#        2.2.3 set environmental variables where rates of change are important in integration ####
tE<-c(tE,list(dMLD  = c(diff(tE$MLD),(tE$MLD[length(tE$MLD)]-tE$MLD[1]))
             ,dCICE = c(diff(tE$CICE),(tE$CICE[length(tE$CICE)]-tE$CICE[1]))
      )) # end concatenation to tE


#    2.3 Timestep vectors invariant over integration - save time here ####

pTaxa<-c("pDi","pSm")
Jmax<-lapply(pTaxa,function(s,a,Temp){
       fnPh_PhotosynthesisMaxJeffery(a[[s]]$Consume$params$MuMax,Temp)
       },a,tE$Temp_MLD)
names(Jmax)<-pTaxa

Ji<-lapply(pTaxa,fnPh_MeanGrowthRateInTstep,a,cE,tE,Jmax) # J with insolation limitation
names(Ji)<-pTaxa


ratio_FeC<-sapply(names(a) # read all iron to carbon ratios for converting carbon to iron
                  ,function(s,a){if(is.null(a[[s]]$Attr$r_FeC)) NA else a[[s]]$Attr$r_FeC},a)
ratio_SiC<-sapply(names(a) # read all silica to carbon ratios for converting carbon to silica
                  ,function(s,a){if(is.null(a[[s]]$Attr$r_SiC)) NA else a[[s]]$Attr$r_SiC},a)

sigmaPB<-sapply(names(a) # read all P/B ratios and estimate sigma for each taxon in estimation vector (not have PB means not used)
                ,function(s,a){if(is.null(a[[s]]$Attr$PBratio)) NA else fnFindSigma(a[[s]]$Attr$PBratio)},a)


# generate parameters for Holling functional response and add to tV
# strip depth ranges from pool data and generate overlap
# strip ingestion rate and selectivity from pool data
# combine availability * selectivtiy = vulnerability

H_Ingest<-list(I = NULL  # maximum ingestion rate matrix - cols[consumer] rows[resource]
               ,v = NULL) # vulnerability matrix cols[consumer] rows[resource] - availability * selectivity

tV <-list(pDi = list(Jmax = Jmax[["pDi"]]
                    ,Ji = Ji[["pDi"]]
                    ) # end pDi list
         ,pSm = list(Jmax = Jmax[["pSm"]]
                    ,Ji = Ji[["pSm"]]
                    ) # end pSm list
         ,nRatio = list(r_FeC = ratio_FeC # note ratio name is the same used in attributes
                       ,r_SiC = ratio_SiC   )
         ,nDissolve = matrix(c(rep(0.0,nStepsPerYr) 
                               ,0.015*exp(0.06*tE$Temp_MLD)^2
                               ,0.015*exp(0.06*tE$Temp_Deep)^2
                               ,rep(0.005,nStepsPerYr)               
                               ,rep(0.0,nStepsPerYr)
                               ,0.015*exp(0.06*tE$Temp_MLD)^2
                               ,0.015*exp(0.06*tE$Temp_Deep)^2
                               ,rep(0.005,nStepsPerYr)  
                               ,rep(0.0,nStepsPerYr)
                               ,pmin(1.32E16*exp(-11200/(tE$Temp_MLD+273.15)),0.02) # Hauck et al 2013
                               ,pmin(1.32E16*exp(-11200/(tE$Temp_Deep+273.15)),0.02) # Hauck et al
                               ,rep(0.005,nStepsPerYr)),nrow=nStepsPerYr,dimnames=list(NULL,  
                         c("dCaSI","dCaM","dCaD","dCaS","dFeSI","dFeM","dFeD","dFeS"
                           ,"dSiSI","dSiM","dSiD","dSiS")))
         ,est = list(mu = sigmaPB^2, sigma=sigmaPB)
        ) # end tV

# set up X0 from input data

X0<-sapply(names(a),function(s,a){a[[s]]$X0},a)

eX<-X0[!is.na(tV$est$mu)]

# call optim here

res<-optim(eX, fnSOweb_estimator, X0=X0,nTimeSteps=nTimeSteps,a=a,cE=cE,tE=tE,tV=tV)
# res<-optim(res$par, fnSOweb_estimator, X0=X0,nTimeSteps=nTimeSteps,a=a,cE=cE,tE=tE,tV=tV)

# then project and plot

# Test
     #   Xest<-X0[!is.na(tV$est$mu)] # for trial projection.  Delete when Xest is estimated
Xest<-res$par
res<-fnSOweb_project(Xest,X0,nTimeSteps,a,cE,tE,tV)

head(res)

plotX<-function(p,t,X,Title){
  Ylim<-c(0,max(X[p,],na.rm=TRUE)*1.1)
  plot(NA,NA,xlim=c(0,370),ylim=Ylim,type="l",xlab="Days",ylab="Abundance",main=Title)
  sapply(seq(1,length(p),1),function(l,p,X,t){
     lines(t,X[p[l],],col=l)
      },p,X,t)
    }
plotX(c("pDi","pSm"),seq(1,(nStepsPerYr+1),1),res,"Phytoplankton")
plotX(c("nFeM","nFeD","nFeSI"),seq(1,(nStepsPerYr+1),1),res,"Nutrient Fe")
plotX(c("nSiM","nSiD","nSiSI"),seq(1,(nStepsPerYr+1),1),res,"Nutrient Si")
plotX(c("dFeM","dFeD","dFeSI"),seq(1,(nStepsPerYr+1),1),res,"Detrital Fe")


# plot(seq(1,(nStepsPerYr+1),1),res["nFeM",],ylim=c(-10,10))
# plot(seq(1,(nStepsPerYr+1),1),res["nFeM",],ylim=c(-10,10),type="l")
# plot(seq(1,(nStepsPerYr+1),1),res["nFeSI",],ylim=c(0,2),type="l")
# plot(seq(1,(nStepsPerYr+1),1),res["pSm",],ylim=c(-10,10),type="l")
# plot(seq(1,(nStepsPerYr+1),1),res["pSm",],ylim=c(-0.5,2),type="l")
# plot(seq(1,(nStepsPerYr+1),1),res["pSm",],ylim=c(-0.5,1.2),type="l")

