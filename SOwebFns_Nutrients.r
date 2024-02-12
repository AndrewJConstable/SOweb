# SOwebFns_Nutrients.r

# by Andrew Constable
# December 2023
# last update 20240206

# 1. Nutrients main functions ####

fnN_Consume<-function(s,sParams,X,a,cE,tE,tV,tStep){ # primary production
  res<-X*0
  vE = list(Dmax = cE$DepthMax
               ,MLD = tE$MLD[tStep]
               ,dMLD = tE$dMLD[tStep]
               ,CICE = tE$CICE[tStep]
               ,dCICE = tE$dCICE[tStep]
               ,Dissolve = tV$nDissolve[tStep,]
               ) # end depths
  dPool<-sapply(sParams$pool,function(pool,Act,s,X,a,cE,vE){
          return(do.call(Act[[pool]]$fn,list(pool,Act[[pool]]$params,s,X,a,cE,vE)))
          },sParams$actions,s,X,a,cE,vE)
  res[sParams$pool]<-as.vector(dPool)
    return(res)
} # end fnN_Consume


fnN_Produce<-function(s,C,sParams,X,a,cE,tE,tV,tStep){ # production based on iron consumption
  return(sum(C[,s])) # vector of amount of each pool consumed in units X
} # end fnN_Produce

# 2. Nutrients supporting functions ####

# changeMLD means the depth pool (MLD or Deep) that increases in depth will gain 
#    conc*dMLD from that pool which is reduced
# i.e. if dMLD is negative then nutrients will be transferred from mixed layer to deep
fnN_changeMLD<-function(pool    # the pool being consumed
                        ,params # parameters for the function of the pool being consumed
                        ,s      # the consumer pool
                        ,X      # vector of abundances (per area)
                        ,a      # attributes of all pools
                        ,cE     # universal constants
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
                       ,a      # attributes of all pools
                       ,cE     # universal constants
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
                            ,a      # attributes of all pools
                            ,cE     # universal constants
                            ,vE     # not used 
){
 return(X[pool]*vE$Dissolve[pool])
} # end fnN_changeMLD

