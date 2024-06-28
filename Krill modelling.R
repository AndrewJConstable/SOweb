# krill modelling

L<-function(a){60*(1-exp(-0.48*a))} # mm
W<-function(L){4E-6*L^3.204} # g
NatAge<-function(a){exp(-0.8*a)}
Nplus<-function(N){N/(1-exp(-0.8))}

a<-c(2:7)
N<-NatAge(a)
N<-N/N[1]
n<-length(a)
N[n]<-Nplus(N[n])
cbind(N,L(a),W(L(a)))
Wmean<-sum(W(L(a))*N)/sum(N) #g
WmnTonne<-Wmean*1E-6


Catch<-8.695E6
KrillNinCatch<-Catch/WmnTonne
KrillNinCatch
