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

  EstWhichX<-!is.na(tV$est$mu)
  
  X <- matrix(NA,length(X0),(nTimeSteps+1),dimnames=list(names(X0),NULL))
  
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

  # divide last values by first values Bend/B0
  Xrel<-X[EstWhichX,(nTimeSteps+1)]/X[EstWhichX,1]
  Mu<-tV$est$mu[EstWhichX]
  Sigma<-tV$est$sigma[EstWhichX]
  
  # use vector list in vE which is mu and sigma for each taxon given PB ratio
  d<-sapply(seq(1,length(Xrel),1),function(x,Xrel,Mu,Sigma){
    return(dlnorm(Xrel,Mu,Sigma)) # probability density for a given X1/X0 i.e. mode=1
    },Xrel,Mu,Sigma)

     return(-sum(log(d))) # return the negative log likelihood
   
      } # end fnSOweb_estimator

fnSOweb_project<-function(eX,X0,nTimeSteps,a,cE,tE,tV){ # estimator is over a one year period
  # some useful blogs
  #  https://www.r-bloggers.com/2019/08/maximum-likelihood-estimation-from-scratch/
  #  https://www.r-bloggers.com/2013/08/fitting-a-model-by-maximum-likelihood/
  
  EstWhichX<-!is.na(tV$est$mu)
  
  X <- matrix(NA,length(X0),(nTimeSteps+1),dimnames=list(names(X0),NULL))
  
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

