# Old code

fnProjPellaTomlinson<-function(y,B,par,catch,Yrs,res){
  if(y>(Yrs-1)){
    return(res)
  } else {
    B<-B+dP(B,par)-fnM(B,par)-catch
    y<-y+1
    res<-rbind(res,data.frame(y=y,B=B))
    return(fnProjPellaTomlinson(y,B,par,catch,Yrs,res))
  }
} # end fnProj


fnTimeToRecover<-function(Bfrom,Bto,par,M,YrsMax=500,Tol=1E-2){  # maximum years for projection
  Brecovered<-Bto*(1-Tol)
  B<-Bfrom
  y<-1
  while (B>0 & B<Brecovered & y<YrsMax) {
    Bp<-B
    B<-B*(1-M)+dP(B,par)
    if(B<0) B<-0
    y<-y+1
  }
  if(B==0 | y==YrsMax) RecoveryTime<-Inf else RecoveryTime<-(Brecovered-Bp)/(B-Bp)+y-1
  return(RecoveryTime)
} # end fnTimeToRecover

# old 

fnProjYear<-function(B,pP,pC,pF=0){
  dP<-dP(B$Bp,pP$par)
  if(pC$useHolling) cC<-B$Bc*pC$QBmax*fnHolling(B$Bp,B$Bc,pC$Holling) else cC<-fnM(B$Bp,pP$par)
  cF<-pF
  print(c(B$Bp,dP,cC,cF,B$Bc,cC*pC$A,B$Bc*(pC$RB+pC$M)))
  Bp1<-B$Bp+dP-cC-cF
  Bc1<-B$Bc+cC*pC$A-B$Bc*(pC$RB+pC$M)
  return(list(Bp=Bp1,Bc=Bc1))
}
Year<-seq(0,50,1)
res<-data.frame(B)
for(i in Year[-1]){
  B<-fnProjYear(B,pP,pC,pF)
  res<-rbind(res,B)
}
res<-cbind(Year,res)


# r for production to match loss at given biomass
findR_loss<-function(r,B,Loss,par){par$r<-r;return((dP(B,par)-Loss)^2)} 


findPhiRB0<-function(par,Ppar,RecoveryYr,ePB,B0propMax){
  #Ppar$phi<-par[1]
  Ppar$r<-par[1]
  Bmax_dP<-nlm(findBmax_dP,0.5*Ppar$K,Ppar)$estimate  # guess if for phi=1
  #B0<-Bmax_dP+(Ppar$K-Bmax_dP)*par[3]
  B0<-nlm(findKrillB0,((Bmax_dP+Ppar$K)/2),Ppar,Bmax_dP,dP(Bmax_dP,Ppar),B0propMax)$estimate
  print(c(Bmax_dP,B0))
  M<-ePB
  Time<-fnTimeToP0(Bmax_dP,(B0*M),Ppar,M)
  print(Time)
  return((Time-RecoveryYr)^2)
} # end findRB0

#Bmax_dP<-nlm(findBmax_dP,0.5*Krill$par$K,Krill$par)$estimate  # guess if for phi=1
#Rmax<-Inf
#Rmin<-nlm(findR_loss,1,Krill$Bmax_dP,Krill$Bmax_dP*Krill$ePB,Krill$par)$estimate # perfect solution if it is possible



findPhi<-function(par,Ppar,ePB,PredKasPropKrillMaxProd,Yrs){
  Ppar$phi<-par
  Bmax_dP<-nlm(findBmax_dP,0.5*Ppar$K,Ppar)$estimate  # guess if for phi=1
  B0<-nlm(findKrillB0,((Bmax_dP+Ppar$K)/2),Ppar,Bmax_dP,dP(Bmax_dP,Ppar),PredKasPropKrillMaxProd)$estimate
  Rmax<-nlm(findR_loss,1,B0,dP(B0,Ppar),Ppar)$estimate # perfect solution if it is possible
  Rmin<- (-log(Bmax_dP*(Ppar$K/B0-1)/(Ppar$K-Bmax_dP))/Yrs) # logistic without loss
  Ppar$r<-optimize(find_R_recovery,interval=c(Rmin,Rmax)
                   ,Ppar=Ppar,Bmax_dP=Bmax_dP,B0=B0,M=dP(B0,Ppar)/B0,RecoveryYr=Yrs)$minimum
  return((dP(B0,Ppar)/B0-ePB)^2)
} # end findPhi





# in order to achieve an ePB of the value here, need to find the shape of the production curve first i.e. phi
# That routine determines the main parameters here.  But they are redetermined here.

#Krill$par$phi<-optimize(findPhi,interval=c(0.1,10),Ppar=Krill$par,ePB=Krill$ePB,PredKasPropKrillMaxProd=Pred$KasPropKrillMaxProd,Yrs=Krill$Recovery$Years)$minimum

# with Phi then solve for the following:



Rmax<-Inf
Rmin<-nlm(findR_loss,1,Krill$Bmax_dP,Krill$Bmax_dP*Krill$ePB,Krill$par)$estimate # perfect solution if it is possible
Krill$par$r<-optimize(find_R_recovery,interval=c(Rmin,Rmax)
                      ,Ppar=Krill$par,Bmax_dP=Krill$Bmax_dP,B0=Krill$B0,M=dP(Krill$B0,Krill$par)/Krill$B0,RecoveryYr=Krill$Recovery$Years)$minimum

##############################################
#######################################################################
# Scenario 1 : no change in Carrying Capacity

plot(NA,NA,xlim=Xlim,ylim=Ylim,xlab="Year",ylab="Population Status", main="Constant Carrying Capacity")
Yrs<-c(Xlim[1]:Xlim[2])

# lines to reflect where the decision rules relate to over time in reference to the initial reference point at time 0

RefPt_neutral<-matrix(c(Xlim,1,1),ncol=2)
Target_neutral<-matrix(c(Xlim,RefPt_neutral[,2]*0.75),ncol=2)
Critical_neutral<-matrix(c(Xlim,RefPt_neutral[,2]*0.2),ncol=2)

lines(RefPt_neutral[,1],RefPt_neutral[,2],lty=4,col="red")
lines(Target_neutral[,1],Target_neutral[,2],lty=2,col="red")
lines(Critical_neutral[,1],Critical_neutral[,2],lty=3,col="red")

# vertical lines separating the different time periods in the simulation

Yrs<-c(Xlim[1]:Xlim[2])
Yr_Fishing_start<-10
Yr_Fishing_end<-70
Yr_recovery_min<-Yr_Fishing_end+Time_for_recovery_min
Yr_recovery_max<-Yr_recovery_min+Time_for_recovery_max-Time_for_recovery_min
abline(v=c(Yr_Fishing_start,Yr_Fishing_end,Yr_recovery_min,Yr_recovery_max),lty=c(3,3,5,5))

# stock trajectory

B<-K_RefPt_0
Yield<-dP(K_target_0*K_RefPt_0,K_RefPt_0,K_r_0,phi)

for (i in 2:length(Yrs)) {
  deltaB<-dP(B[(i-1)],K_RefPt_0,K_r_0,phi)
  if(Yrs[i]>=Yr_Fishing_start && Yrs[i]<=Yr_Fishing_end) Catch<-Yield else Catch<-0
  #  print(c(Yrs[i],deltaB,Catch))
  
  B<-c(B,(B[(i-1)]+deltaB-Catch)) 
}


lines(Yrs,(B/K_RefPt_0),lty=1)



#######################################################################
# Scenario 2 : declining carrying capacity with no change in yield (target levels and critical levels are adjusted as proportion)

plot(NA,NA,xlim=Xlim,ylim=Ylim,xlab="Year",ylab="Population Status", main="Declining Carrying Capacity")
Yrs<-c(Xlim[1]:Xlim[2])


# lines to reflect where the decision rules relate to over time in reference to the initial reference point at time 0
Kgrad<- (-0.01)
Yr_startDecline<-40
RefPt<-rep(K_RefPt_0,length(Yrs))
RefPt[Yrs>Yr_startDecline]<-K_RefPt_0+Kgrad*(Yrs[Yrs>Yr_startDecline]-Yr_startDecline)
Target<-RefPt*K_target_0
Critical<-RefPt*K_critical_0

lines(Yrs,RefPt,lty=4,col="red")
lines(Target,lty=2,col="red")
lines(Critical,lty=3,col="red")

# vertical lines separating the different time periods in the simulation

Yrs<-c(Xlim[1]:Xlim[2])
Yr_Fishing_start<-10
Yr_Fishing_end<-70
Yr_recovery_min<-Yr_Fishing_end+Time_for_recovery_min
Yr_recovery_max<-Yr_recovery_min+Time_for_recovery_max-Time_for_recovery_min
abline(v=c(Yr_Fishing_start,Yr_Fishing_end,Yr_recovery_min,Yr_recovery_max),lty=c(3,3,5,5))

# stock trajectory

B<-K_RefPt_0
Yield<-dP(K_target_0*K_RefPt_0,K_RefPt_0,K_r_0,phi)

for (i in 2:length(Yrs)) {
  deltaB<-dP(B[(i-1)],RefPt[i],K_r_0,phi)
  if(Yrs[i]>=Yr_Fishing_start && Yrs[i]<=Yr_Fishing_end) Catch<-Yield else Catch<-0
  #  print(c(Yrs[i],deltaB,Catch))
  
  B<-c(B,(B[(i-1)]+deltaB-Catch)) 
}


lines(Yrs,(B/B[1]),lty=1)

#######################################################################
# Scenario 3 : rapid shift in carrying capacity and then stable

plot(NA,NA,xlim=Xlim,ylim=Ylim,xlab="Year",ylab="Population Status", main="Rapid shift in Carrying Capacity")
Yrs<-c(Xlim[1]:Xlim[2])


# lines to reflect where the decision rules relate to over time in reference to the initial reference point at time 0
Kgrad<- (-0.05)
Yr_startDecline<-40
Yr_endDecline<-45
RefPt<-rep(K_RefPt_0,length(Yrs))
RefPt[Yrs>Yr_startDecline]<-K_RefPt_0+Kgrad*(Yrs[Yrs>Yr_startDecline]-Yr_startDecline)
RefPt[Yrs>Yr_endDecline]<-RefPt[sum(Yrs<=Yr_endDecline)]

Target<-RefPt*K_target_0
Critical<-RefPt*K_critical_0


lines(Yrs,RefPt,lty=4,col="red")
lines(Target,lty=2,col="red")
lines(Critical,lty=3,col="red")

# vertical lines separating the different time periods in the simulation

Yrs<-c(Xlim[1]:Xlim[2])
Yr_Fishing_start<-10
Yr_Fishing_end<-70
Yr_recovery_min<-Yr_Fishing_end+Time_for_recovery_min
Yr_recovery_max<-Yr_recovery_min+Time_for_recovery_max-Time_for_recovery_min
abline(v=c(Yr_Fishing_start,Yr_Fishing_end,Yr_recovery_min,Yr_recovery_max),lty=c(3,3,5,5))

# stock trajectory

B<-K_RefPt_0
Yield<-dP(K_target_0*K_RefPt_0,K_RefPt_0,K_r_0,phi)

for (i in 2:length(Yrs)) {
  deltaB<-dP(B[(i-1)],RefPt[i],K_r_0,phi)
  if(Yrs[i]>=Yr_Fishing_start && Yrs[i]<=Yr_Fishing_end) Catch<-Yield else Catch<-0
  B<-c(B,(B[(i-1)]+deltaB-Catch)) 
}


lines(Yrs,(B/B[1]),lty=1)

#######################################################################
# Scenario 4 : catstrophic shift of stock abundance due to change in food web - leading to change in critical level

plot(NA,NA,xlim=Xlim,ylim=Ylim,xlab="Year",ylab="Population Status", main="Rapid shift in Carrying Capacity & Catastrophic change")
Yrs<-c(Xlim[1]:Xlim[2])


# lines to reflect where the decision rules relate to over time in reference to the initial reference point at time 0
Kgrad<- (-0.05)
Yr_startDecline<-40
Yr_endDecline<-45
RefPt<-rep(K_RefPt_0,length(Yrs))
RefPt[Yrs>Yr_startDecline]<-K_RefPt_0+Kgrad*(Yrs[Yrs>Yr_startDecline]-Yr_startDecline)
RefPt[Yrs>Yr_endDecline]<-RefPt[sum(Yrs<=Yr_endDecline)]

Target<-RefPt*K_target_0
Critical<-RefPt*K_critical_0

# catastrophic shift in stock status at 50 years and change in critical level
CritGrad<- (0.0005)
Yr_catastrophe<-50
Critical[Yrs>Yr_catastrophe]<-1.2*Critical[sum(Yrs<=Yr_catastrophe)]+CritGrad*(Yrs[Yrs>Yr_catastrophe]-Yr_catastrophe)



lines(Yrs,RefPt,lty=4,col="red")
lines(Target,lty=2,col="red")
lines(Critical,lty=3,col="red")

# vertical lines separating the different time periods in the simulation

Yrs<-c(Xlim[1]:Xlim[2])
Yr_Fishing_start<-10
Yr_Fishing_end<-70
Yr_recovery_min<-Yr_Fishing_end+Time_for_recovery_min
Yr_recovery_max<-Yr_recovery_min+Time_for_recovery_max-Time_for_recovery_min
abline(v=c(Yr_Fishing_start,Yr_Fishing_end,Yr_recovery_min,Yr_recovery_max),lty=c(3,3,5,5))

# stock trajectory

B<-K_RefPt_0
Yield<-dP(K_target_0*K_RefPt_0,K_RefPt_0,K_r_0,phi)

for (i in 2:length(Yrs)) {
  deltaB<-dP(B[(i-1)],RefPt[i],K_r_0,phi)
  if(Yrs[i]>=Yr_Fishing_start && Yrs[i]<=Yr_Fishing_end) Catch<-Yield else Catch<-0
  if(Yrs[i]==Yr_catastrophe) deltaB<-(-0.5*B[length(B)])
  B<-c(B,(B[(i-1)]+deltaB-Catch)) 
}


lines(Yrs,(B/B[1]),lty=1)

