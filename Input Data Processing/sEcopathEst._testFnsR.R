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
  { Iwt<-I/QB
    if(Iwt<1.5) Iwt<-1 else Iwt<-(Iwt/1.5)^0.5
    res<-Iwt*sum( (I*res/(I+sumRes) - d*QB)^2)}
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

# Test functions ####

root<-"/Users/acon/Desktop/_w/_p/SOfoodweb & EcoVelocity/Ecosystem network synthesis/Ecopath files/Rinputs/"
rootEcopathEst<-"/Users/acon/Desktop/_w/_r/SOweb/EcopathEst/"
# Inputs for Ecopath Analysis ####

units<-365 # time units: 1=year, 365=days

iUseLog<-FALSE # log transform in root function

#Ecopath_Mods<-c("Ballerini2014","Ballerini2014_Stand","Dahood2019","Gurney2014_Yr2000","Gurney2014_Stand","Hill2021","Hill2021_Stand","Maldonado2016_Yr1900","Maldonado2016_Yr2008","McCormack2019","Pinkerton2010","Pinkerton2010_Stand","Subramanium2020","Surma2014")

Ecopath_Mods<-c("Pinkerton2010")


EcoMod<-Ecopath_Mods[1]
# designate Ecopath model in order to generalise code ####

eMod<-fnReadEcopath(root,EcoMod,ReplaceNA=TRUE)

# diet - columns are the proportion diet composition for a consumer
# check diet adds to 1 for each predator (not always the case) ####

In_B<-eMod$params[,"B"]
In_PB<-eMod$params[,"P.B"]
In_QB<-eMod$params[,"Q.B"]/units
In_EE<-eMod$params[,"EE"]
Diet<-eMod$diet


res<-fnFindPool_Iv(13,In_B,In_QB,Diet,iUseLog)









# Check how well Ecopath model is consistent in diet matrix ####

fnCheckEcoPathFit(In_B,In_QB,In_PB,PredDiet$d,In_EE)

# estimate vulnerability matrix (v) and predator max ingestion rates (I) ####


#     test on a single predator ####
tgt<-17
units<-365 # time units: 1=year, 365=days
mask<-Diet[,tgt]>0
iB<-In_B[mask]
iQB<-In_QB[tgt]/units
id<-Diet[mask,tgt]
#iUseLog<-FALSE # log transform in root function

Vstart<-c(iQB/2,rep(0.9,length(id)))
boundLo<-1E-10  
boundHi<-c(Inf,rep(1,length(id)))  # does not like Inf in vector but if it was to be changed to 1 then need to ensure that estimation of I is between 0 and 1, which may require scaling
rootFnIhat(Vstart,B=iB,QB=iQB,d=id,useLog=iUseLog)

Vstart<-c(iQB/2,rep(0.2,length(id)))

res<-optim_ARS(par=c(iQB/2,rep(0.5,length(id))),fn=rootFnIhat,lower=boundLo,upper=boundHi
               ,iter=10000,iter_adapt=10,trace=FALSE,B=iB,QB=iQB,d=id,useLog=iUseLog)


res<-trialV(0.1,iQB,iB,iQB,id,iUseLog,boundLo,boundHi)
res<-searchV(iB,iQB,id,iUseLog,boundLo,boundHi)
res1<-minV(res)

# test result ####
tgt_j<-17 # predator to target
EcoMod<-Ecopath_Mods[1]
eMod<-fnReadEcopath(root,Ecopath_Mods[1],ReplaceNA=TRUE)
load(paste0(rootEcopathEst,"EcopathEst_",EcoMod,"_day_OLS.Rdata"))
tB<-eMod$params[,"B"]
tI<-EcopathEst$IngestionMax[tgt_j]
tv<-tv<-EcopathEst$v_matrix[,tgt_j]
tgt_i<-which(tv>0)

print(tv)
xMax<-max(tB)*1.5
tx<-seq(0,xMax,xMax*0.01)
fnIprop<-function(i,I,v,B){v[i]*B[i]/(I+sum(v*B))}
plot(tx,(tx/(tI+tx)),type="l")
points(tB[tgt_i],fnIprop(tgt_i,tI,tv,tB))
tv_alt<-tv
tv_alt[tB>2]<-0
points(tB[tgt_i],fnIprop(tgt_i,tI,tv_alt,tB),col="red")


print(tB[tgt_i])



