# SOwebFn_Detritus.r

# by Andrew Constable
# December 2023
# last update 20240206

# 1. Detritus main functions

fnD_Consume<-function(s,sParams,X,a,cE,tE,tV,tStep){ # primary production
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
} # end fnD_Consume

fnD_Produce<-function(s,C,sParams,X,a,cE,tE,tV,tStep){ 
# note: production can be moving a type of detritus from one pool to another
  # requiring only a simple transfer from consumption to production.
  # iron and silica are 'produced' from the accumulation of carbon detritus 
  # from living pools or carcasses.  The ratio of iron or silica to carbon in the living/carcass pools
  # is then used to produce the iron or silica.
  
#  i. same units - mole to mole or carbon to carbon
#  ii. carbon mass to mole of consumer
  tV<-c(tV,list(dCarbonPools_C = C[,c("dCaSI","dCaM","dCaD","dCaS")]))
  
  res<-sapply(sParams$pool,function(pool,vC,Act,s,a,cE,tV){
    return(do.call(Act[[pool]]$fn,list(pool,Act[[pool]]$params,vC,s,a,cE,tV)))
  },C[,s],sParams$actions,s,a,cE,tV)
  return(sum(as.vector(res)))
    } # end fnD_Produce

# 2. Detritus Supporting functions

fnD_MLtoDeep<-function(pool    # the pool being consumed
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
  if(vE$dMLD<0) d<- (-vE$dMLD)/vE$MLD else d<-0
  return(X[pool]*(d+params$rSink/vE$MLD))
} # end fnN_changeMLD

fnD_DeepToML<-function(pool    # the pool being consumed
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
  if(vE$dMLD>0) return(X[pool]/(vE$Dmax-vE$MLD)*vE$dMLD) else return(0)
} # end fnN_changeMLD


# changeSI - same approach as for nutrients 
fnD_changeSI<-function(pool    # the pool being consumed
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



fnD_consumePhMortality<-function(pool,Params,s,X,a,cE,vE){
  return(X[pool]*Params$Ragg)
} # end fnD_consumePhMortality
  
fnD_consumeFeScavenge<-function(pool,Params,s,X,a,cE,vE){
  return(min(X[Params$CarbonScavengingPool]*cE$Fe_scavenge$Rate*X[pool],X[pool]*cE$Fe_scavenge$FreePropDFe))
}

fnD_sinkFromDeep<-function(pool,Params,s,X,a,cE,vE){
  return(Params$rSink/(cE$DepthMax-vE$MLD)*X[pool])
} # end fnD_sinkFromDeep

fnD_resuspend<-function(pool,Params,s,X,a,cE,vE){
  return(Params$rSuspend*X[pool])
} # end fnD_resuspend

fnD_produceSameUnits<-function(pool,Params,vC,s,a,cE,tV){ # return quantity
  return(vC[pool]) 
} # end fnD_produceSameUnits

fnD_produceMoleFromCarbonPool<-function(pool,Params,vC,s,a,cE,tV){ #  pool is carbon
  return(tV$dCarbonPools_C[pool,Params$Cpool]*cE$C_MassToMole*tV$nRatio[[a[[s]]$Attr$Which_C_ratio]][[pool]]) # convert carbon to mole based on ratio of pool
  } # end fnD_produceMoleFromCarbon

