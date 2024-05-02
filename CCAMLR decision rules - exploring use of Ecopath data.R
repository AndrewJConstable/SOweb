library(ggplot2)
K<-0.8
cP<-0.8E-6
cM<-0.1
cR<-0.05

Krill<-list(Holling = list(p50 = 0.5, q = 0))
Whale<-list(QBhat=1,Holling = list(p50 = 0.25, q = 1,c=1,g=3,doPC=TRUE))

Bc<-seq(0.1,1,0.1)
Kc<-seq(0,1.5,0.1)
Hg<-seq(1,3,0.25)

# examine effect of Hg
res<-do.call(rbind,lapply(Hg,function(g,K,B,H){H$g<-g
                         res<-fnHolling(K,B,H)
                         res<-as.data.frame(cbind(rep(g,length(B)),B,res))
                         names(res)<-c("g","B","Crate")
                         return(res)
                         },Kc[length(Kc)],Bc,Whale$Holling))
res$g<-factor(res$g)
p<-ggplot(res,aes(x=B,y=Crate,colour=g))+geom_line()
p

# examine effect of Bc given Hg
g<-4
res<-do.call(rbind,lapply(Bc,function(B,g,K,H){H$g<-g
res<-fnHolling(K,B,H)
res<-as.data.frame(cbind(rep(B,length(K)),K,res))
names(res)<-c("B","K","Crate")
return(res)
},g,Kc,Whale$Holling))
res$B<-factor(res$B)
p<-ggplot(res,aes(x=K,y=Crate,colour=B),title=paste0("G = ",g))+geom_line()
p

PredCompetition<-function(Bp,Bc,H){if(H$doPC) return(1/(Bp/(Bc*H$c))^H$g) else return(0)}

fnHolling<-function(Bp,Bc,H){Bp^(H$q+1)/(H$p50^(H$q+1)+(Bp*(1+PredCompetition(Bp,Bc,H)))^(H$q+1))}




fnScaleQB<-function(Bc,Kp,QBhat,c){
  x<-Bc*QBhat/Kp
  x[x>1]<-1
  xR<-(x-c$x50)
  adjY<-(1-c$y50)/0.5
  xR[xR>0]<-xR[xR>0]/adjY
  
  r<-1/(1+exp(-c$K*xR))

  r[xR<=0]<-r[xR<=0] * c$y50/0.5
  r[xR>0] <- 1-(1-r[xR>0]) * adjY  
  return(r)}

Bc<-seq(0,3,0.01)
Kp<-1.5
plot(Bc,fnScaleQB(Bc,Kp,Whale$QBhat,Whale$QBscaleParams),type="l")

fnQ.B<-function(B,K,Qhat,phi){r<-B*0
                              useB<-B/K<1
                              r[useB]<-Qhat*B[useB]*(1-B[useB]/K)^phi
                              return(r)}

dBc<-function(Bc,Bp,QBhat,U,RBhat,Mhat,omega,Hp50,Hq){
              (QBhat*fnHolling(Bp,Hp50,Hq)*(1-U)-RBhat-Mhat)*Bc
              }# end fn

B<-seq(0,1,0.01)
plot(seq(0,1,0.01),fnQ.fnII(seq(0,1,0.01),0.01),type="l")
plot(seq(0,1,0.01),fnQ.B(seq(0,1,0.01),1,0.2,0.2),type="l")
#plot(B,fnP(B 
,B*fnQ.B(B,K,cQ*fnQ.fnII(0.1,0.05),0.5),0.8,cR,cM),type="l")

cQ<-0.28
plot(B,fnP(B,B*cQ*fnQ.fnII(0.1,0.05),0.8,cR,cM),type="l")
