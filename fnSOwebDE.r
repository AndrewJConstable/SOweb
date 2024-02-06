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

# 1. Inputs ####

#    1.1. Minimisation inputs ####
# time0 relates to 1st element of a temporal vector.  The matrix of state variables 
nYrs<-1
nStepsPerYr<-365 
day<-c(1:nStepsPerYr)

Date1<-as.Date("2015-07-01")


#    1.2. Input Parameters ####

a<-list(  Nu = NULL
         ,Ph = list(Di = NULL
                   ,Sm = NULL)
         ,Zo = NULL
         ,De = NULL
       )


#    1.3. Input Environmental Data ####

fPathIn<-"./IO in/"
fPathOut<-"./IO out/"
load(paste0(fPathIn,"Fe_means.Rdata")) # df Iron concentration rows(MEASO area) cols(n, name, value)
load(paste0(fPathIn,"Si_means.Rdata")) # df Silicate concentration rows(MEASO area) cols(n, name, value)
load(paste0(fPathIn,"Lat_means.Rdata")) # df Latitude mean for MEASO areas rows(MEASO area) cols(n, name, value)
load(paste0(fPathIn,"Insol_means.Rdata"))  # matrix Insolation W.m-2.d-1 rows(months) cols(MEASO area)
load(paste0(fPathIn,"MLD_means.Rdata"))  # matrix Mixed Layer Depth m rows(months) cols(MEASO area)
load(paste0(fPathIn,"CICE_means.Rdata"))  # matrix Sea ice percent cover % rows(months) cols(MEASO area) 
load(paste0(fPathIn,"Temp_MLD_means.Rdata"))  # matrix Temperature mixed layer oC rows(months) cols(MEASO area) 
load(paste0(fPathIn,"Temp_Deep_means.Rdata"))  # matrix Temperature MLD to 1000m oC rows(months) cols(MEASO area) 

TargetArea <-"AOA"

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
nTimeSteps<-nYrs*nStepsPerYr
DaysInMonth<-c(31,28,31,30,31,30,31,31,30,31,30,31)
DayMidMonth<-c(15,14,15,15,15,15,15,15,15,15,15,15)
#    2.2 Environmental variables ####
#        2.2.1 changing order of environmental variables to start with the first month ####
Month1<-format(Date1,"%m")

fnReorderMonthVector<-function(m,V){
       e<-which(names(eTemp_MLD)%in%paste0("M",m))
       if(e==1) return(V) else return(V[c(e:12,1:(e-1))]) 
       } # end fn

eInsol     <-fnReorderMonthVector(Month1,eInsol)
eMLD       <-fnReorderMonthVector(Month1,eMLD)
eCICE      <-fnReorderMonthVector(Month1,eCICE)
eTemp_MLD  <-fnReorderMonthVector(Month1,eTemp_MLD)
eTemp_Deep <-fnReorderMonthVector(Month1,eTemp_Deep)

#        2.2.2 interpolation of environmental vectors to time steps ####
# interpolation to all time steps based on the monthly values being mid month
# do one year and then replicate for the number of years needed

# vector for days corresponding to months plus beginning and end of year for interpolation
DayMM<-c(0,DayMidMonth[1],(cumsum(DaysInMonth)[1:11]+DayMidMonth[2:12]),366)

fnInterpolateVectorToTimeSteps<-function(V,VtSteps,tSteps){
if((length(VtSteps)-length(V))==2)  {V_0<-(V[1]+V[length(V)])/2; V<-c(V_0,V,V_0)}
approx(VtSteps,V,tSteps)$y
}

eInsol_tSteps     <-rep(fnInterpolateVectorToTimeSteps(eInsol,DayMM,day),nYrs)
eMLD_tSteps       <-rep(fnInterpolateVectorToTimeSteps(eMLD,DayMM,day),nYrs)
eCICE_tSteps      <-rep(fnInterpolateVectorToTimeSteps(eCICE,DayMM,day),nYrs)
eTemp_MLD_tSteps  <-rep(fnInterpolateVectorToTimeSteps(eTemp_MLD,DayMM,day),nYrs)
eTemp_Deep_tSteps <-rep(fnInterpolateVectorToTimeSteps(eTemp_Deep,DayMM,day),nYrs)

#        2.2.3 set environmental variables where rates of change are important in integration ####
dMLD_tSteps <- c(diff(eMLD_tSteps),(eMLD_tSteps[length(eMLD_tSteps)]-eMLD_tSteps[1]))
dCICE_tSteps <- c(diff(eCICE_tSteps),(eCICE_tSteps[length(eCICE_tSteps)]-eCICE_tSteps[1]))



#    2.3 Timestep vectors invariant over integration - save time here ####

tIV <-list(    Nu = NULL
              ,Ph = list(Di = list(Jmax = )
                         ,Sm = NULL)
              ,Zo = NULL
              ,De = NULL
)

Jmax
Jbar


## RHS of NPZD system (from Melbourne-Thomas et al 2015) ####

fnSOwebDE <- function(k # vector element to read  (not used by JMT)
              ,x # X vector - NPZD
              ,u # mortality rate
              ,a # parameters
              ,mld # mixed layer depth for time step
              ,no3e # nitrate concentration for time step
              ,w # change in mixed layer depth for time step
              ,J # phytoplankton growth rate
){ # start function
  ## Grazing
  G <- (a[9]*a[10]*x[2]^2)/(a[9] + a[10]*x[2]^2)
  wp <- max(w,0)
  c(## Nitrate equation
    a[13]*x[4]+a[12]*x[3]-J*x[2]+(wp*(no3e-x[1]))/mld,
    ## Phytoplankton equation
    (J-wp/mld)*x[2]-a[6]*x[2]-a[7]*x[2]^2-G*x[3],
    ## Zooplankton equation
    (a[8]*G-a[12]-w/mld)*x[3]-a[11]*x[3]^2,
    ## Detritus equation
    (1-a[8])*G*x[3]+a[11]*x[3]^2+a[6]*x[2]+a[7]*x[2]^2-(a[13]+(wp+a[14])/mld)*x[4])
}

