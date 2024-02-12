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
  dPool<-sapply(sParams$pool,function(pool,Act,s,X,cE,vE){
    return(do.call(Act[[pool]]$fn,list(pool,Act[[pool]]$params,s,X,a,cE,vE)))
  },sParams$actions,s,X,a,cE,vE)
  res[sParams[[s]]$pool]<-as.vector(dPool)
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

fnD_consumePhMortality<-function(pool,Params,s,X,a,cE,vE){
  return(X[pool]*Params$Ragg)
} # end fnD_consumePhMortality
  
fnD_consumeFeScavenge<-function(pool,Params,s,X,a,cE,vE){
  return(min(X[Params$CarbonScavengingPool]*cE$Fe_scavenge$Rate*X[pool],X[pool]*cE$Fe_scavenge$FreePropDFe))
}

fnD_sinkFromML<-function(pool,Params,s,X,a,cE,vE){
  return(Params$rSink/vE$MLD*X[pool])
  } # end fnD_sinkFromML

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

