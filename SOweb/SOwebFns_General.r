# SOwebFn_General.r

# by Andrew Constable
# December 2023
# last update 20240206


# estimator ####

# P/B ratio is standard deviation of lognormal distribution with mode=1 
fnFindSigma<-function(SD){
  fHolling<-function(a,x){x^a[1]/(a[2]^a[1]+x^a[1])}
  guess<-fHolling(c(1.683762,1.637206),SD) # from exploration of relationship
  fnF<-function(v,SD){return((SD-exp(1.5*v)*sqrt(exp(v)-1))^2)}
  v<-nlm(fnF,guess,SD)
  return(sqrt(v$estimate))}


fnSOweb_estimator<-function(eX,X0,nTimeSteps,a,cE,tE,tV){ # estimator is over a one year period
  # some useful blogs
#  https://www.r-bloggers.com/2019/08/maximum-likelihood-estimation-from-scratch/
#  https://www.r-bloggers.com/2013/08/fitting-a-model-by-maximum-likelihood/
print(eX)
  # if any of eX are negative then return Inf
  
if(sum(eX<=0)>0){ res<-Inf} else { # all positive
  EstWhichX<-!is.na(tV$est$mu)
  
  X <- matrix(NA,length(X0),(nTimeSteps+1),dimnames=list(names(X0),NULL))

# modify X0 nutrients and detritus pools based on cE$Conc_t0  
  whichL<-which(cE$Conc_t0_layer=="M")
  X0[names(cE$Conc_t0_layer)[whichL]]<-unlist( cE$Conc_t0[names(cE$Conc_t0_layer)[whichL]] )*tE$MLD[1]
  whichL<-which(cE$Conc_t0_layer=="D")
  X0[names(cE$Conc_t0_layer)[whichL]]<-unlist( cE$Conc_t0[names(cE$Conc_t0_layer)[whichL]] )*(cE$DepthMax-tE$MLD[1])
  whichL<-which(cE$Conc_t0_layer=="SI")
  X0[names(cE$Conc_t0_layer)[whichL]]<-unlist( cE$Conc_t0[names(cE$Conc_t0_layer)[whichL]] ) # assuming sea ice is 1m thick
  
  X[,1]<-X0
  X[EstWhichX,1]<-eX
  for(k in 2:(nTimeSteps+1)) { # 2:n is based on the matrix of state variables
    X[,k]<-X[,k-1]+fnSOwebDE((k-1) # vector element to read  (not used by JMT)
                             ,X[,(k-1)] # X vector
                             ,a   # parameters
                             ,cE  # environmental constants
                             ,tE  # time step vectors of environmental variables
                             ,tV) # time step vectors for use as needed in integration
  } # end do loop

  # divide last values by first values Bend/B0
  Xrel<-X[EstWhichX,(nTimeSteps+1)]/X[EstWhichX,1]
  Mu<-tV$est$mu[EstWhichX]
  Sigma<-tV$est$sigma[EstWhichX]
  
  # use vector list in vE which is mu and sigma for each taxon given PB ratio
  d<-sapply(seq(1,length(Xrel),1),function(x,Xrel,Mu,Sigma){
    return(dlnorm(Xrel,Mu,Sigma)) # probability density for a given X1/X0 i.e. mode=1
    },Xrel,Mu,Sigma)
  res<--sum(log(d))
} # if all eX positive

     return(res) # return the negative log likelihood
   
      } # end fnSOweb_estimator

fnSOweb_project<-function(eX,X0,nTimeSteps,a,cE,tE,tV){ # estimator is over a one year period
  # some useful blogs
  #  https://www.r-bloggers.com/2019/08/maximum-likelihood-estimation-from-scratch/
  #  https://www.r-bloggers.com/2013/08/fitting-a-model-by-maximum-likelihood/
  
  EstWhichX<-!is.na(tV$est$mu)
  
  X <- matrix(NA,length(X0),(nTimeSteps+1),dimnames=list(names(X0),NULL))

    # modify X0 nutrients and detritus pools based on cE$Conc_t0  
  whichL<-which(cE$Conc_t0_layer=="M")
  X0[names(cE$Conc_t0_layer)[whichL]]<-unlist( cE$Conc_t0[names(cE$Conc_t0_layer)[whichL]] )*tE$MLD[1]
  whichL<-which(cE$Conc_t0_layer=="D")
  X0[names(cE$Conc_t0_layer)[whichL]]<-unlist( cE$Conc_t0[names(cE$Conc_t0_layer)[whichL]] )*(cE$DepthMax-tE$MLD[1])
  whichL<-which(cE$Conc_t0_layer=="SI")
  X0[names(cE$Conc_t0_layer)[whichL]]<-unlist( cE$Conc_t0[names(cE$Conc_t0_layer)[whichL]] ) # assuming sea ice is 1m thick
  
  
  X[,1]<-X0
  X[EstWhichX,1]<-eX
  for(k in 2:(nTimeSteps+1)) { # 2:n is based on the matrix of state variables
    cat("Day ",(k-1),"\n",sep="")
    
    X[,k]<-X[,k-1]+fnSOwebDE((k-1) # vector element to read  (not used by JMT)
                             ,X[,(k-1)] # X vector
                             ,a   # parameters
                             ,cE  # environmental constants
                             ,tE  # time step vectors of environmental variables
                             ,tV) # time step vectors for use as needed in integration
  } # end do loop
  
return(X) # return the negative log likelihood
  
} # end fnSOweb_project


fnSOwebDE <- function(t # vector element to read  (not used by JMT)
                      ,X # X vector - NPZD
                      ,a # parameters
                      ,cE
                      ,tE # new - time step vectors of environmental variables
                      ,tV # new - time step vectors for use as needed in integration
){ # start function
  # Generate X by X matrix of consumption - rows(consumed) cols(consumer) 
  # usual consumption of predators and prey
  # include consumption of detritus by nutrients
  Action<-"Consume"
  Consumption<-sapply(names(X),function(s,X,a,cE,tE,tV,t){
    if(is.null(a[[s]][[Action]])) return(0) else if(a[[s]][[Action]]$fn=="" | is.null(a[[s]][[Action]]$fn)) return(X*0) else
      return(do.call(a[[s]][[Action]]$fn,list(s,a[[s]][[Action]]$params,X,a,cE,tE,tV,t))) # return vector of consumed taxa
  },X,a,cE,tE,tV,t)
  Consumption<-as.matrix(Consumption)
  
# crude check and correction for over consumption of a pool ####
#   if iron is overconsumed then reduction in iron also should be applied to silica #
#   if silica is overconsumed then reduction in silica must be applied to silica consuming phytoplankton
#   for higher trophic levels - need to think this through

  #???????????? hardwired for the moment 
  checkC<-apply(Consumption,1,sum,na.rm=TRUE)
  XoverC<-checkC>X
  count<-0
  while(sum(XoverC)>0){
    count<-count+1
    pools<-names(which(XoverC)) # pools that have been overconsumed
    if(sum(pools=="nSiM")>0){ # adjust silicate first as need to reduce primary production 
        Adj<-X["nSiM"]/checkC["nSiM"]
        Consumption["nFeM","pDi"]<-Consumption["nFeM","pDi"]*Adj
        Consumption["nSiM",]<-Consumption["nSiM",]*Adj
    } else # Silicate check 
      if(sum(pools=="nFeM")>0){ # adjust iron second and adjust silicate consumption in diatoms 
        Adj<-X["nFeM"]/checkC["nFeM"]
        Consumption["nSiM","pDi"]<-Consumption["nSiM","pDi"]*Adj
        Consumption["nFeM",]<-Consumption["nFeM",]*Adj
      } else {
        for(i in c(1:length(pools))) Consumption[pools[i],]<-Consumption[pools[i],]*X[pools[i]]/checkC[pools[i]]
        }
    checkC<-apply(Consumption,1,sum,na.rm=TRUE)
    XoverC<-checkC>X
  }  

    # Vector of Production from consumption
  
  Action<-"Produce"
  Production<-sapply(names(X),function(s,X,a,cE,tE,tV,k,C){
    if(is.null(a[[s]][[Action]])) return(0) else if(a[[s]][[Action]]$fn=="" | is.null(a[[s]][[Action]]$fn)) return(0) else
      return(do.call(a[[s]][[Action]]$fn,list(s,C,a[[s]][[Action]]$params,X,a,cE,tE,tV,k))) # return vector of consumed taxa
  },X,a,cE,tE,tV,k,Consumption)
  Production<-as.vector(Production)
  
  # Update X
  return(Production-apply(Consumption,1,sum,na.rm=TRUE))
} # end DE

fnInterpolateVectorToTimeSteps<-function(V,VtSteps,tSteps){
  if((length(VtSteps)-length(V))==2)  {V_0<-(V[1]+V[length(V)])/2; V<-c(V_0,V,V_0)}
  approx(VtSteps,V,tSteps)$y
}

