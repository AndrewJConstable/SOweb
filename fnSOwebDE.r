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

# 1. Minimisation inputs

TargetArea <-"AOA"

# 2. Input functions and data ####


source("SOwebFns.r")
source("SOwebData.r")

# 3. Subset Environmental Data ####

nFe_0      <-Fe_means[Fe_means[,"MEASO_area"]==TargetArea,3]
nSi_0      <-Si_means[Si_means[,"MEASO_area"]==TargetArea,3]
eLat       <-Lat_means[Lat_means[,"MEASO_area"]==TargetArea,3]
eInsol     <-Insol_means[,TargetArea]
eMLD       <-MLD_means[,TargetArea]
eCICE      <-CICE_means[,TargetArea]
eTemp_MLD  <-Temp_MLD_means[,TargetArea]
eTemp_Deep <-Temp_Deep_means[,TargetArea]


# 2. Pre-minimisation simplification and standardisation ####

#    2.1 Managing time ####
day<-c(1:nStepsPerYr)
nTimeSteps<-nYrs*nStepsPerYr

#    2.2 Environmental variables ####
#
#        2.2.1 constants

cE$kw <- 1/cE$kw  # light attenuation in water - modified for computational efficiency

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
                  ,function(s,a){a[[s]]$Attr$r_FeC},a) 
names(ratio_FeC)<-names(a)

tV <-list(pDi = list(Jmax = Jmax[["pDi"]]
                    ,Ji = Ji[["pDi"]]
                    ) # end pDi list
         ,pSm = list(Jmax = Jmax[["pSm"]]
                    ,Ji = Ji[["pSm"]]
                    ) # end pSm list
         ,nRatio = list(r_FeC = ratio_FeC # note ratio name is the same used in attributes
                       ,r_SiC = NULL   )
        ) # end tIV


# set up X

X<-rep(1000,length(a))
names(X)<-names(a)

## RHS of NPZD system (from Melbourne-Thomas et al 2015) ####

fnSOwebDE <- function(k # vector element to read  (not used by JMT)
              ,X # X vector - NPZD
              ,a # parameters
              ,cE
              ,tE # new - time step vectors of environmental variables
              ,tV # new - time step vectors for use as needed in integration
#              ,mld # mixed layer depth for time step
#              ,no3e # nitrate concentration for time step
#              ,w # change in mixed layer depth for time step
#              ,J # phytoplankton growth rate
              ){ # start function
Xp1<-X*0


# Generate X by X matrix of consumption - rows(consumed) cols(consumer) 
   # usual consumption of predators and prey
   # include consumption of detritus by nutrients
      Action<-"Consume"
      Consumption<-sapply(names(X),function(taxon,X,a,cE,tE,tV){
         if(a[[taxon]][[Action]]$fn=="" | is.null(a[[taxon]][[Action]]$fn) | is.null(a[[taxon]][[Action]])) return(X*0) else
          return(do.call(a[[taxon]]$Consume$fn,list(taxon,a[[taxon]]$Consume$params),X,a,cE,tE,tV,k)) # return vector of consumed taxa
      },X,a,cE,tE,tV)
      Consumption<-as.matrix(do.call(cbind,Consumption))
      dimnames(Consumption)<-list(names(X),names(X))

# Vector of Mortality from consumption

# Vector of Growth from consumption

# Generate matrix of detrital accumulation (natural mortality) - rows(taxa) cols(taxa)
#    note only detrital pools will have values > 0

# Generate vector of import - X 

# Generate vector of export - X

Nutrients<-1 # just for testing

Xp1[a$Ph$Di$X] <- X[a$Ph$Di$X]*(tV$Ph$Di$Ji[k]*Nutrients)
Xp1[a$Ph$Sm$X] <- X[a$Ph$Sm$X]*(tV$Ph$Sm$Ji[k]*Nutrients)

    ## Grazing
#       G <- (a[9]*a[10]*x[2]^2)/(a[9] + a[10]*x[2]^2)
#       wp <- max(w,0)
  return(Xp1)
   #c(## Nitrate equation
         # a[13]*x[4]+a[12]*x[3]-J*x[2]+(wp*(no3e-x[1]))/mld,
    ## Phytoplankton equation
         # (J-wp/mld)*x[2]-a[6]*x[2]-a[7]*x[2]^2-G*x[3],
    ## Zooplankton equation
         # (a[8]*G-a[12]-w/mld)*x[3]-a[11]*x[3]^2,
    ## Detritus equation
         # (1-a[8])*G*x[3]+a[11]*x[3]^2+a[6]*x[2]+a[7]*x[2]^2-(a[13]+(wp+a[14])/mld)*x[4]
} # end DE


X <- matrix(NA,2,(nTimeSteps+1))

X[,1]<-c(100,50)

for(k in 2:(nTimeSteps+1)) { # 2:n is based on the matrix of state variables
#  JA <- Jmax[k-1]*X[1,k-1]/(a[5]+X[1,k-1])
#  JB <- min(J[k-1],JA)
#  X[,k] <- X[,k-1] + dt*F(k-1,X[,k-1],u[k-1],a,mld[k-1],no3e[k-1],w[k-1],JB)
X[,k]<-X[,k-1]+fnSOwebDE((k-1) # vector element to read  (not used by JMT)
                        ,X[,(k-1)] # X vector - NPZD
                        ,a # parameters
                        ,tE # new - time step vectors of environmental variables
                        ,tV) # new - time step vectors for use as needed in integration
  }

tx<-seq(0,(nTimeSteps/365),1/365)
plot(tx,log10(X[1,]),type="l",col="blue")#,ylim=c(0,10000))
lines(tx,log10(X[2,]),col="red")

