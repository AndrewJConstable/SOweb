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

fnH_consumeHolling2_KY<-function(pool,Params,s,X,a,cE,vE){ # after Koen-Alonso & Yodzis, 2005
   # s is the consumer, pool is the resource
  return((vE$H_Ingest$I[s,pool]*vE$H_Ingest$v[s,pool]*X[pool])/ # numerator
          (vE$H_Ingest$I[s,pool]+sum(vE$H_Ingest$v[,s]*X)))     # denominator
} # end fnH_consumeHolling2_KY

