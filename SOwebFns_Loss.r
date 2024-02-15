# SOwebFn_Loss.r
# consumption based on loss functions across taxa e.g. respiration, other mortality
# maybe included in specific pools later

# by Andrew Constable
# December 2023
# last update 20240206

# 1. Loss main functions

fnL_Consume<-function(s,sParams,X,a,cE,tE,tV,tStep){ # primary production
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
} # end fnL_Consume

fnL_Produce<-function(s,C,sParams,X,a,cE,tE,tV,tStep){ 
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

# 2. Loss Supporting functions

