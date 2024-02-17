# generate size ane mdf per cohort if needed
for(s in c(1:length(Pool_Forage))){
if(!is.null(Pool_Forage[[s]]$cohort$fn)){
  if(is.null(Pool_Forage[[s]]$cohort$size)) Pool_Forage[[s]]$cohort$size<-Pool_Forage[[s]]$cohort$fn$size(Pool_Forage[[s]]$cohort$params)
  if(is.null(Pool_Forage[[s]]$cohort$freq)) Pool_Forage[[s]]$cohort$freq<-Pool_Forage[[s]]$cohort$fn$freq(Pool_Forage[[s]]$cohort$params)
  } # end estimate size per cohort
} # end do loop on size per cohort



# estimate size and weight-freq of cohorts given 
# Linf, age at 95% Linf, Weight-Length params (a,b), AgeMax

fnCohortSummary<-function(A0=0,AgeMax=10,A_1pc=10,A_95Linf=5,Linf=100,t0=(-0.1),WL=c(1,3),L95){
  M <- (-log(0.01)/A_1pc)
  K=log(0.05)/(-A_95Linf+t0)
  a<-seq(A0,AgeMax,1)
  N<-exp(-M*a)
  N[length(a)]<-N[length(a)]/(1-exp(-M))
  N<-N/sum(N)
  L<-fnLength_vB(Linf,K,t0,a)
  W<-WL[1]*L^WL[2]
  Freq<-N*W/sum(N*W)
data.frame(a,N,L,W,Freq)  
} # end fnCohortSummary

fnLength_vB<-function(Linf,K,t0,a){Linf*(1-exp(-K*(a-t0)))}
  

tmp<-fnCohortSummary()

