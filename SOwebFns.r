# SOwebFns.r

# by Andrew Constable
# December 2023
# last update 20240206


# General ####

fnInterpolateVectorToTimeSteps<-function(V,VtSteps,tSteps){
  if((length(VtSteps)-length(V))==2)  {V_0<-(V[1]+V[length(V)])/2; V<-c(V_0,V,V_0)}
  approx(VtSteps,V,tSteps)$y
}

fnPh_MeanGrowthRateInTstep<-function(s,a,cE,tE,Jmax){
  # mean growth rate is the weighted mean (by percent coverage of sea ice) of the daily mean growth rate in sea ice and in open water 
  
  # open water
  JepOpen<-EvansParslow(s,a,cE,tE,cE$par$w,Jmax)  
  
  # sea ice
  JepSI<-EvansParslow(s,a,cE,tE,cE$par$si,Jmax)  
  
  # return weighted mean
  return(JepOpen*tE$OpenWater+JepSI*tE$CICE)
} # end function

# Phytoplankton ####

fnPh_PhotosynthesisMaxJeffery<-function(MuMax,T){MuMax*exp(0.06*T)}

#     Evans-Parslow routine modified from Melbourne-Thomas et al 2015
#          Angle of incident radiation is dependent on time of day and day of the year (tilt of the Earth) and influences the 
#          incident radiation reaching depths.  Corrections for this are included in the description on page 66 in Eq 5.
 
EvansParslow <- function( # JMT Equation 5 Average Growth Rate given latitude, day length and mixed layer depth
                          # note Jmax(T) from Jeffery et al
   s     # character name for phytoplankton group in list
  ,a      # biological parameters
  ,cE     # environmental constants
  ,tE     # environmental time series
  ,PAR    # while PAR is present as a list if cE, only one is parsed.  The one parsed is either with or without sea ice present.
  ,Jmax   # list of timestep vectors of Jmax for different phytoplankton groups  
                 ) { # start function
#  using the following data:
#   tE$day    # Julian day of the year
#  ,tE$Insol    # solar radiation at surface (daily forcing variable)
#  ,tE$MLD    # mixed layer depth (daily forcing variable)
#  ,a$Ph[[Ph]]$alpha  # photosynthetic efficiency
#  ,cE$kw     # constant (Table 1) - light attenuation in water   
#  ,Jmax[[Ph]]   # maximum growth rate as 'nutrient limited (light saturated) growth'
#  ,cE$Lat    # latitude
  ##    
  ## Compute daily averaged light-limited growth rate 
  ## analytical integration over layer thickness and day after Evans and
  ## Parslow (1985)
  #fx1 <- 2*pi*(day+192)/365
  fx1 <- 2*pi*tE$day/365
  declin <- 0.006918-  
    (0.399912*cos(fx1)+0.006758*cos(2*fx1)+0.002697*cos(3*fx1)) +
      +(0.070257*sin(fx1)+0.000907*sin(2*fx1)+0.001480*sin(3*fx1))

  ## Compute solar angle at noon (and assume that this is the
  ## equivalent daily averaged incidence angle for direct+diffuse
  ## radiation). cobeta is cos(incidence angle of solar radiation at
  ## noon)
  fx1 <- pi/180*cE$Lat
  cobeta <- pmax(sin(fx1)*sin(declin)+cos(fx1)*cos(declin),0)
  cobeta <- sqrt(1-(1-cobeta^2)/1.33^2)

  ## Length of day
  fx2 <- pmax(pmin(-tan(fx1)*tan(declin),1),-1)
  daylen <- pmax(1.0e-12,acos(fx2)/pi)

  rayb <- pmax(1.0e-15,exp(-tE$MLD/(cobeta*cE$kw)))

  daylen <- daylen/2
  radbio <- pmax(1.0,PAR*tE$Insol)
  vpbio <- Jmax[[s]]
  fx1 <- daylen^2*vpbio/(a[[s]]$Consume$params$alpha*radbio)
  fx3 <- fx1
  fx4 <- fx1/rayb
  fu1 <- sqrt(fx3^2+daylen^2)
  fu1 <- fu1-daylen*log((daylen+fu1)/fx3)
  fu2 <- sqrt(fx4^2+daylen^2)
  fu2 <- fu2-daylen*log((daylen+fu2)/fx4)
  J <- -2*vpbio/tE$MLD*cobeta*cE$kw*(fu1-fu2-fx3+fx4)  # note cE$kw is 1/attenuation coefficient

  J
}

fnPhConsume<-function(s,sParams,X,a,cE,tE,tV,tStep){ # primary production
  res<-X*0
  MLD<-tE$MLD[tStep]
  PhyConc<-X[s]*cE$C_MassToMole/MLD # convert mass to C mole conc
  poolConc<-sapply(sParams$pool,function(f,X,MLD){X[f]/MLD},X,MLD)
  names(poolConc)<-sParams$pool
      # determine nutrient uptake rates, take minimums & adjust J
  J<-tV[[s]]$Ji[tStep]*min(sapply(sParams$pool,function(f,FC,P){FC[f]/(FC[f]+P$k[P$pool%in%f])},poolConc,sParams),na.rm=TRUE)
  J[is.na(J)]<-0
    res[sParams$pool]<-sapply(sParams$pool,function(f,s,J,a,tV,PhyMoleC){
                       return(J*PhyMoleC*tV$nRatio[[a[[f]]$Attr$Which_C_ratio]][[s]])
                       },s,J,a,tV,PhyConc*MLD)
  return(res) # vector of amount of each pool consumed in units X
} # end fnPhConsume


fnPhProduceFe<-function(s,C,sParams,X,a,cE,tE,tV,tStep){ # production based on iron consumption
  C_FeRatio<-1/a[[s]]$Attr$r_FeC 
  return(C["nFeM",s]*cE$C_MoleToMass*C_FeRatio) # vector of amount of each pool consumed in units X
} # end fnPhProduceFe


# Nutrients ####

fnN_Consume<-function(s,sParams,X,a,cE,tE,tV,tStep){ # primary production
  res<-X*0
  vE = list(Dmax = cE$DepthMax
               ,MLD = tE$MLD[tStep]
               ,dMLD = tE$dMLD[tStep]
               ,CICE = tE$CICE[tStep]
               ,dCICE = tE$dCICE[tStep]
               ,Dissolve = tV$nDissolve[tStep,]
               ) # end depths
  dPool<-sapply(sParams$pool,function(pool,Act,s,X,vE){
          return(do.call(Act[[pool]]$fn,list(pool,Act[[pool]]$params,s,X,vE)))
          },sParams$actions,s,X,vE)
  print(dPool)
  res[sParams[[s]]$pool]<-as.vector(dPool)
    return(res)
} # end fnN_Consume


fnN_Produce<-function(s,C,sParams,X,a,cE,tE,tV,tStep){ # production based on iron consumption
  return(sum(C[,s])) # vector of amount of each pool consumed in units X
} # end fnN_Produce


# changeMLD means the depth pool (MLD or Deep) that increases in depth will gain 
#    conc*dMLD from that pool which is reduced
# i.e. if dMLD is negative then nutrients will be transferred from mixed layer to deep
fnN_changeMLD<-function(pool    # the pool being consumed
                        ,params # parameters for the function of the pool being consumed
                        ,s      # the consumer pool
                        ,X      # vector of abundances (per area)
                        ,vE # list of depth and sea ice values for function
                                   # ,Dmax   # maximum water depth
                                   # ,MLD    # mixed layer depth in time step
                                   # ,dMLD   # change in mixed layer depth in time step
                        ){
  if(vE$dMLD>0 & params$do_dMLDgt0) {
    return(X[pool]/(vE$Dmax-vE$MLD)*vE$dMLD)
   } else if(vE$dMLD<0 & !params$do_dMLDgt0) {
     return(X[pool]/vE$MLD*(-vE$dMLD))
   } else return(0) 
} # end fnN_changeMLD
  
# changeSI - same approach as for MLD
fnN_changeSI<-function(pool    # the pool being consumed
                        ,params # parameters for the function of the pool being consumed
                        ,s      # the consumer pool
                        ,X      # vector of abundances (per area)
                        ,vE # list of depth and sea ice values for function
                                    # ,CICE   # proportion of area covered by sea ice
                                    # ,dCICE  # change in proportion of area covered by sea ice
){
    if(vE$dCICE>0 & params$do_dCICEgt0) {
    return(X[pool]/vE$MLD*vE$dCICE/(1-vE$CICE)) # assume freezing of only top 1 m
  } else if(vE$dCICE<0 & !params$do_dCICEgt0) {
        return(X[pool]*(-vE$dCICE)/vE$CICE)
  } else return(0) 
} # end fnN_changeMLD

fnN_remineraliseD<-function(pool    # the pool being consumed
                            ,params # parameters for the function of the pool being consumed
                            ,s      # the consumer pool
                            ,X      # vector of abundances (per area)
                            ,vE     # not used 
){
 return(X[pool]*vE$Dissolve[pool])
} # end fnN_changeMLD

# Detritus ####


fnD_Consume<-function(s,sParams,X,a,cE,tE,tV,tStep){ # primary production
  res<-X*0
  vE = list(Dmax = cE$DepthMax
            ,MLD = tE$MLD[tStep]
            ,dMLD = tE$dMLD[tStep]
            ,CICE = tE$CICE[tStep]
            ,dCICE = tE$dCICE[tStep]
            ,Dissolve = tV$nDissolve[tStep,]
  ) # end depths
  dPool<-sapply(sParams$pool,function(pool,Act,s,X,vE){
    return(do.call(Act[[pool]]$fn,list(pool,Act[[pool]]$params,s,X,vE)))
  },sParams$actions,s,X,vE)
  res[sParams[[s]]$pool]<-as.vector(dPool)
  return(res)
} # end fnD_Consume

fnD_Produce<-function(s,C,sParams,X,a,cE,tE,tV,tStep){ 
# note: three types of production
#  i. same units - mole to mole or carbon to carbon
#  ii. carbon mass to mole of consumer

  res<-sapply(sParams$pool,function(pool,vC,Act,s,a,cE,tV){
    return(do.call(Act[[pool]]$fn,list(pool,Act[[pool]]$params,vC,s,a,cE,tV)))
  },C[,s],sParams$actions,s,cE,tV)
  return(sum(as.vector(res)))
    } # end fnD_Produce

fnD_produceSameUnits<-function(pool,Params,vC,s,a,cE,tV){ # return quantity
  return(vC[pool]) 
} # end fnD_produceSameUnits

fnD_produceMoleFromCarbon<-function(pool,Params,vC,s,a,cE,tV){ #  pool is carbon
  return(vC[pool]*cE$C_MassToMole*tV$nRatio[[a[[s]]$Attr$Which_C_ratio]][[pool]]) # convert carbon to mole based on ratio of pool
  } # end fnD_produceMoleFromCarbon

