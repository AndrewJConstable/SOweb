# SOwebFn_General.r

# by Andrew Constable
# December 2023
# last update 20240206

## RHS of NPZD system (from Melbourne-Thomas et al 2015) ####

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

