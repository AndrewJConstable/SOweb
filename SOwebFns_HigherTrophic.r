# SOwebFns_HigherTrophic.r

# by Andrew Constable
# December 2023
# last update 20240206

# 1. Higher Trophic pools main functions ####
  
fnH_Consume<-function(s,sParams,X,a,cE,tE,tV,tStep){ # primary production
  res<-X*0
  # calc sum vHat*Bi
  
  Frate<-sapply(sParams$pools,sParams$fn etc.....)
  
  res[sParams$pool]<-Frate*X[sParams$pool]
  return(res) # vector of amount of each pool consumed in units X
} # end fnH_Consume


fnH_Produce<-function(s,C,sParams,X,a,cE,tE,tV,tStep){ 

  return() # production of X[s]
} # end fnH_Produce

# 2. Higher Trophic pools supporting functions ####

fnH_consumeHolling2<-function(pool,Params,s,X,a,cE,vE){

    return(res)
}

