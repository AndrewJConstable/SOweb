# Script for processing Ecopath files

# Libraries to support optimisation ####

library(PopED)

# Functions ####

fnReadEcopath<-function(root,fname,ReplaceNA=TRUE){
  names<-read.csv(paste(root,fname,"_names.csv",sep=""),header=TRUE)
  params<-read.csv(paste(root,fname,"_params.csv",sep=""),header=TRUE)
  diet<-read.csv(paste(root,fname,"_diet.csv",sep=""),header=FALSE)
  
  if(ReplaceNA){
    params[is.na(params[,"B"]),"B"]<-0
    params[is.na(params[,"P.B"]),"P.B"]<-1 # consume biomass even though no production occurs (e.g. detritus)
    params[is.na(params[,"Q.B"]),"Q.B"]<-0
    diet[is.na(diet)]<-0
  }
  return(list(names=names,params=params,diet=diet))
}

calcPredConsume<-function( # returns vector of prey consumed
  I   # scalar or vector - max ingestion rates for each prey of predator
  ,Br  # vector - biomasses of prey
  ,v   # vector - vulnerability of prey to predator
  ,Bp  # scalar - Biomass of predator
){
  return(sum(Holling2(Br,I,v)*I)*Bp)
}

Holling2<-function(   # returns vector proportion of max ingestion consumed
  Br # vector - resource (prey) biomasses
  ,I  # scalar
  ,v  # vector - vulnerability of prey to predator
){
  mask<-v>0
  res<-v*0
  sum_vBr<-sum(v[mask]*Br[mask])
  res[mask]<-v[mask]*Br[mask]/(I+sum_vBr)
  return(res)
}

fnCheckEcoPathFit<-function(B,QB,PB,d,EE){
  n<-length(B)
  pool<-seq(1,n,1)
  res<-lapply(pool,function(i,B,QB,PB,d,EE){
    left<-sum(unlist(
      sapply(pool,function(j,i,d,B,QB)
      {d[i,j]*B[j]*QB[j]},i,d,B,QB)
    )) # end unlist sum
    right<-EE[i]*PB[i]*B[i]
    return((left-right)/right)
  },B,QB,PB,d,EE)
  return(unlist(res))
} # end check function


#     Root function for minimisation routine ####

rootFnIhat<-function(par,B,QB,d,useLog){
  I<-par[1]
  v<-par[-1]
  res<-v*B
  sumRes<-sum(res)
  if(useLog) {res<-sum( (log(I*res/(I+sumRes)) - log(d*QB))^2)} else
  {res<-sum( (I*res/(I+sumRes) - d*QB)^2)}
  return(res)
} # end find Ihat

#     Other functions to support minimising I,v ####

trialV<-function(v,I,B,QB,d,iUseLog,boundLo,boundHi){
  maxRun<-400
  resI<-NULL
  resV=NULL
  resValue=NULL
  # unadjusted
  Vstart<-c(I,rep(v,length(d)))
  res<-optim_ARS(par=Vstart,fn=rootFnIhat,lower=boundLo,upper=boundHi,iter=10000,iter_adapt=10,max_run=maxRun,trace=FALSE,B=B,QB=QB,d=d,useLog=iUseLog)
  resI<-c(resI,res$par[1])
  resV<-rbind(resV,res$par[-1])
  resValue<-c(resValue,res$ofv)
  
  # unadjusted 0.3*I
  Vstart<-c((I*0.3),rep(v,length(d)))
  res<-optim_ARS(par=Vstart,fn=rootFnIhat,lower=boundLo,upper=boundHi,iter=10000,iter_adapt=10,max_run=maxRun,trace=FALSE,B=B,QB=QB,d=d,useLog=iUseLog)
  resI<-c(resI,res$par[1])
  resV<-rbind(resV,res$par[-1])
  resValue<-c(resValue,res$ofv)
  
  testI<-resI[which(resValue==min(resValue,na.rm=TRUE))]
  
  # adjust Iup
  Vstart<-c(testI*1.5,resV[1,])
  res<-optim_ARS(par=Vstart,fn=rootFnIhat,lower=boundLo,upper=boundHi,iter=10000,iter_adapt=10,max_run=maxRun,trace=FALSE,B=B,QB=QB,d=d,useLog=iUseLog)
  resI<-c(resI,res$par[1])
  resV<-rbind(resV,res$par[-1])
  resValue<-c(resValue,res$ofv)
  
  # adjust Idown
  Vstart<-c(testI*0.5,resV[1,])
  res<-optim_ARS(par=Vstart,fn=rootFnIhat,lower=boundLo,upper=boundHi,iter=10000,iter_adapt=10,max_run=maxRun,trace=FALSE,B=B,QB=QB,d=d,useLog=iUseLog)
  resI<-c(resI,res$par[1])
  resV<-rbind(resV,res$par[-1])
  resValue<-c(resValue,res$ofv)
  
  # adjust inverse V
  Vstart<-c(testI,(1-resV[1,]))
  res<-optim_ARS(par=Vstart,fn=rootFnIhat,lower=boundLo,upper=boundHi,iter=10000,iter_adapt=10,max_run=maxRun,trace=FALSE,B=B,QB=QB,d=d,useLog=iUseLog)
  resI<-c(resI,res$par[1])
  resV<-rbind(resV,res$par[-1])
  resValue<-c(resValue,res$ofv)
  
  # adjust invert last V
  Vstart<-c(res$par[1]*0.5,(1-res$par[-1]))
  res<-optim_ARS(par=Vstart,fn=rootFnIhat,lower=boundLo,upper=boundHi,iter=10000,iter_adapt=10,max_run=maxRun,trace=FALSE,B=B,QB=QB,d=d,useLog=iUseLog)
  resI<-c(resI,res$par[1])
  resV<-rbind(resV,res$par[-1])
  resValue<-c(resValue,res$ofv)
  
  return(list(name=c("unadj","unadj_third_I","Iup","Idown","inverseV","lastInvert"),I=resI,v=resV,value=resValue))  
} # end function

minV<-function(r){
  minR<-which(r$value==min(r$value))
  return(list(name=r$name[minR],I = r$I[minR],v = r$v[minR,],value = r$value[minR]))
} # return minimum result from list of results

searchV<-function(B,QB,d,UseLog,boundLo,boundHi){
  res<-lapply(seq(0.1,0.9,0.1),trialV,QB,B,QB,d,UseLog,boundLo,boundHi)
  resMin<-unlist(lapply(res,function(r){
    r1<-which(unlist(r$value)==min(unlist(r$value) ))
    return(r1[1])}))
  resI<-unlist(sapply(seq(1,length(resMin),1),function(i,r,rMin){r[[i]]$I[rMin[i]]},res,resMin))
  resV<-do.call(rbind,lapply(seq(1,length(resMin),1),function(i,r,rMin){r[[i]]$v[rMin[i],]},res,resMin))
  resValue<-unlist(sapply(seq(1,length(resMin),1),function(i,r,rMin){r[[i]]$value[rMin[i]]},res,resMin))
  resName<-unlist(sapply(seq(1,length(resMin),1),function(i,r,rMin){r[[i]]$name[rMin[i]]},res,resMin))
  return(list(name=resName,I=resI,v=resV,value=resValue))  
} # end function

fnFindPool_Iv<-function(j,iB,iQB,Diet,iUseLog){
  print(paste0(EcoMod," Pool ",j))
  mask<-Diet[,j]>0
  nPrey<-sum(mask)
  B  <- iB[mask]
  QB <- iQB[j]
  d  <- Diet[mask,j]
  if(nPrey==0 | sum(B>0)==0) return(c(QB,Diet[,j]))
  boundLo <- c(1E-10,rep(1E-10,length(d)))
  boundHi <- c(Inf,rep(1,length(d)))
  r1<-minV(searchV(B,QB,d,iUseLog,boundLo,boundHi))
  print(r1)
  if(length(r1$I)>1){
    v<-Diet[,j]*0
    v[mask]<-rep(1,length(d))
    rI<-fnCalcI(which(mask)[1]        # target prey species
                ,iB      # vector of biomasses of all taxa
                ,v       # vulnerabilities of biomasses to predator
                ,QB     # consumption per biomass of predator
    )
    r2<-c(nPrey,r1$value[1],rI,v)  
    print("Minimum of Search returned more than one result.  Vulnerabilities returned were:")
    print(r1$v)
    
  } else {
    v<-Diet[,j]*0
    v[mask]<-r1$v
    r2<-c(nPrey,r1$value,r1$I,v)
  }
  return(r2)
} # end findPool_Iv

fnCheckEstI_V<-function(I              # estimated I vector
                        ,v             # estimated v matrix
                        ,d             # diet matrix
                        ,B             # vector of biomasses
                        ,QB            # vector - consumption by each pool
){
  n<-length(QB)
  pool<-seq(1,n,1)
  # vector of total consumption by predators
  dprime<-lapply(pool,function(j,I,v,d,B,QB){
    if(I[j]==0) return(v[,j]*0)
    res<-Holling2(B,I[j],v[,j])*I[j]
    res<-round(((res-d[,j]*QB[j])/(d[,j]*QB[j]))*100,0)
    res[is.nan(res)]<-0
    return(res)
  },I,v,d,B,QB)
  dprime<-do.call(cbind,dprime)
  
  return(dprime)
} # end check function

fnCalcI<-function(i        # target prey species
                  ,B      # vector of biomasses of all taxa
                  ,v       # vulnerabilities of biomasses to predator
                  ,QBj     # consumption per biomass of predator
){sum(v*B)/(v[i]*B[i]/QBj-1)}

# Run Estimation ####

root<-"/Users/acon/Desktop/_w/_p/SOfoodweb & EcoVelocity/Ecosystem network synthesis/Ecopath files/Rinputs/"

Ecopath_Mods<-c("Pinkerton2010")

# Ecopath_Mods<-c("Ballerini2014","Ballerini2014_Stand","Dahood2019","Gurney2014_Yr2000"
#                ,"Gurney2014_Stand","Hill2021","Hill2021_Stand","Maldonado2016_Yr1900"
#                ,"Maldonado2016_Yr2008","McCormack2019","Pinkerton2010","Pinkerton2010_Stand"
#                ,"Subramanium2020","Surma2014")


units<-365
iUseLog=FALSE

resEcopath<-lapply(Ecopath_Mods,function(EcoMod){
  print(EcoMod)
  eMod<-fnReadEcopath(root,EcoMod,ReplaceNA=TRUE)
  checkPredTotalDiet<-apply(eMod$diet,2,sum,na.rm=TRUE)
  Npool<-nrow(eMod$params)
  
  In_B<-eMod$params[,"B"]
  In_PB<-eMod$params[,"P.B"]
  In_QB<-eMod$params[,"Q.B"]/units
  In_EE<-eMod$params[,"EE"]
  Diet<-eMod$diet
  
  checkEcopathFit<-fnCheckEcoPathFit(In_B,In_QB,In_PB,Diet,In_EE)
  rSeq<-seq(1,Npool,1)
  res<-do.call(cbind,lapply(rSeq,fnFindPool_Iv
                            ,In_B,In_QB,Diet,iUseLog))
  nPrey<-res[1,]
  rootValue <- res[2,]
  estI<-res[3,]
  estV<-res[-c(1:3),]
  m0pd<-(1-unlist(eMod$params[,"EE"]))*unlist(eMod$params[,"P.B"])/units
  M0pd<--log(1-m0pd)
  RNA<-(unlist(eMod$params[,"Q.B"])-unlist(eMod$params[,"P.B"]))/unlist(eMod$params[,"Q.B"])
  
  EcopathEst<-list(
     PoolN              = Npool
    ,QB                 = In_QB
    ,checkPredTotalDiet = checkPredTotalDiet
    ,checkEcopathFit    = checkEcopathFit
    ,IngestionMax       = estI
    ,v_matrix           = estV
    ,nPrey              = nPrey
    ,rootValue          = rootValue
    ,NonPredMortRate    = M0pd
    ,Resp_NonAssim      = RNA
  )
  save(EcopathEst,file=paste0("EcopathEst_",EcoMod,"_day_OLS.Rdata"))
  return(EcopathEst)  
})


