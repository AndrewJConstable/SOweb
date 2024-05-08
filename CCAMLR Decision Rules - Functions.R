# CCAMLR Decision Rules - Functions

# 1.0 General ####


# 2.0 Krill ####

#     2.1 Pella-Tomlinson 1969 ####

dP<-function(B,par){B*par$r/par$phi*(1-(B/par$K)^par$phi)}
  
fnPellaTomlinson_dBdt<-function(t,y,parms){
  catch<-if(length(parms$catch)>1) parms$catch[floor(t)+1] else parms$catch 
  list(dP(y,parms$PT$par)-fnM(y,parms$PT)-catch)
}

fnPlotPellaTomlinson<-function(PT,HollingExamples=NULL,Krange=c(0,1),Ylim=NULL,plotMrate=TRUE,plotMtypeH=TRUE,catch=0){
# create dataframe ####
#    X-axis ####
       B<-seq(Krange[1],Krange[2],1E-2)*PT$par$K 

#    Production without mortality and initialise dataframe & max production vector ####
         LineName<-"Productivity"
         Change  <- dP(B,PT$par)
         Line    <- rep(LineName,length(B))
         pd      <- data.frame(Line,B,Change)

#    Mortality rate
         if(plotMrate){
           LineName<-"M_constant"
           Change  <- PT$par$M*B
           Line    <- rep(LineName,length(B))
           pd      <- rbind(pd,data.frame(Line,B,Change))
           }

#    Mortality - functional relationships
         if(plotMtypeH & !is.null(HollingExamples)){
          pd<-rbind(pd,do.call(rbind,lapply(seq(1,length(HollingExamples),1),function(h,H,B,PT){
            PT$par$Holling<-H[[h]]

            PT$par$useMrate<-FALSE
            PB0<-PT$B0*PT$par$M
            PT$par$Mmax<-nlm(findMmax,PB0,PT,PT$B0,PB0)$estimate
            Change  <- fnM(B,PT)
            Line    <- rep(paste0("M_",names(H)[h]),length(B))
            return(data.frame(Line,B,Change))
          },HollingExamples,B,PT)))           
         } # end if
  pd$Line<-factor(pd$Line,levels=unique(pd$Line))

  Maxima<-pd$B[tapply(pd$Change,pd$Line,function(pdchange){which(pdchange==max(pdchange))})]
  names(Maxima)<-levels(pd$Line)
  Maxima<-Maxima[names(Maxima)!="M_constant"]
           # setup graph as needed ####
    Xlim<-Krange*PT$par$K
    if(is.null(Ylim)) Ylim<-c(0,ceiling(max(pd$Change)/5)*5)
    
    
  # plot graph ####
  p<-ggplot(pd)+labs(Title="Annual Production Curve",x="Biomass",y="Production")
  p<-p+theme(panel.background = element_rect(fill="white",colour = "black")
            ,panel.grid.major = element_blank()
            ,panel.grid.minor = element_blank()
            ,axis.line=element_line(colour="black")
            ,axis.text = element_text(size=14)
            ,axis.title=element_text(size=18)
            )# end theme
  p<- p+scale_x_continuous(breaks=seq(Xlim[1],Xlim[2],length.out=6)
                           ,expand=expansion(mult = c(0,0.05)))
  p<- p+scale_y_continuous(breaks=seq(Ylim[1],Ylim[2],length.out=5),expand=expansion(mult = c(0,0.05)))
  p<- p+guides(
    x = guide_axis(minor.ticks = TRUE),
    y = guide_axis(minor.ticks = TRUE)
  )
  
#  p<- p+geom_vline(xintercept=PT$Bmax_dP,colour="black",linetype="dotted") # bionmass with max production
  p<- p+geom_vline(xintercept=PT$B0,colour="black",linetype="twodash") # B0
  
     p<-p+geom_line(aes(x=B,y=Change,colour=Line))
     p<- p+theme(legend.key=element_blank()) # removes boxes from around lines in legend
     
# now do colours and line types as chosen                   
#  ,colour="darkgreen",linetype="dashed")+geom_vline(xintercept=BatPmaxWithMfn,colour="darkgreen",linetype="dotted")
    p<-p+geom_vline(xintercept=Maxima,colour=c("black","darkblue","darkgreen"),linetype=c("dotted","dotted","dotted"))
return(p)
} # end fn


#     2.2 Mortality function for population modelling ####
fnM<-function(B,PT){if(PT$par$useMrate) return(PT$par$M*B) else 
       return(PT$par$Mmax*(B^PT$par$Holling$q)/((PT$par$Holling$p50)^PT$par$Holling$q+B^PT$par$Holling$q))}

#     2.3 Setup routines to calculate parameters from input data ####
findBmax_dP<-function(b,par){return(-log(dP(b,par)))}  # = 0.5K when phi=1
findKrillB0<-function(par,Ppar,bmax,pmax,p){prod<-dP(par,Ppar); return((prod-p*pmax)^2)}
findProdMaxGivenM<-function(b,PT) {return(-log(dP(b,PT$par)-fnM(b,PT)))}
findMmax<-function(Mmax,PT,B0,pB0){PT$par$Mmax<-Mmax;return((fnM(B0,PT)-pB0)^2)}

# 3.0 Predators ####

PredCompetition <- function(Bp,Bc,H){if(H$doPC) return(1/(Bp/(Bc*H$c))^H$g) else return(0)}
fnHolling       <- function(Bp,Bc,H){Bp^(H$q+1)/(H$p50^(H$q+1)+(Bp*(1+PredCompetition(Bp,Bc,H)))^(H$q+1))}

# scaling predator to prey
# The relationship of predator to prey at equilibrium is dependent on the biomass of krill at equilibrium,
#       its productivity at the biomass (relative to maximum production) and the 
#       consumption rate of the predators relative to its maximum consumption rate

findPredQBmax   <- function(im,qb,Cb0,Pb0,H){return((im*fnHolling(Pb0,Cb0,H)-qb)^2)} # QBmax

# 4.0 Food web ####

#     4.1 Change in food web - dBdt ####

fnFoodWeb_dBdt<-function(t,y,parms){
  pE<-parms$pE
  pP<-parms$pP
  pC<-parms$pC
  catch<-if(length(parms$catch)>1) parms$catch[floor(t)+1] else parms$catch 
  
  # do prey
  if(pP$modelK$changingK) pP$par$K<-pE[(floor(t)+1)]
  PdP<-dP(y[1],pP$par)

  # do predation
  r1<-do.call(rbind,lapply(seq(1,length(pC),1),function(p,pC,y,pP){   
    if(pC[[p]]$useHolling) cC<-y[(1+p)]*pC[[p]]$QBmax*fnHolling(y[1],y[(1+p)],pC[[p]]$Holling) else cC<-fnM(y[1],pP)
    return(data.frame( cC=cC,dB=cC*pC[[p]]$A-y[(1+p)]*(pC[[p]]$RB+pC[[p]]$M) ))
    },pC,y,pP))
  r2<-c((PdP-sum(r1[,"cC"])-catch),as.vector(r1[,"dB"]) )
  names(r2)<-names(y)
  return(list(r2)) # end return list
}


# 5.0 Fishery ####

#     4.1 Find gamma for population assessment ####

find_G_Pop<-function(par,PT,B0,Btarget,projYrs){
  g<-par
  res<-ode(y=B0,times=seq(0,projYrs,1),func=fnPellaTomlinson_dBdt,parms=list(PT=PT,catch=g*B0),method="rk4")
  nYrsGTE0<-sum(!is.nan(res[,2]))-1
  if(nYrsGTE0<projYrs) {
    t<-(-res[(nYrsGTE0+1),2]/(res[(nYrsGTE0+1),2]-res[nYrsGTE0,2]))
    return((-Btarget*(2-(nYrsGTE0+t)/projYrs))^2)} else {
      return((res[nrow(res),2]-Btarget)^2)}
}

#     4.2 Find gamma for food web assessment ####

find_G_FW<-function(par,parms,fwB0,Btarget,projYrs){
  g<-par
  parms$catch<-g*fwB0[1]
  res<-ode(y=fwB0,times=seq(0,projYrs,1),func=fnFoodWeb_dBdt,parms=parms,method="rk4")
  # the first time at Btarget relative to projYrs is to be minimised
  y<-res[,1]
  B<-res[,2]
  nYrsGTE0<-max(y[which(!is.nan(B))])
  
  if(nYrsGTE0==projYrs){  # projection to projYrs successful
    whichBmin<-which(B==min(B))
    if(length(whichBmin)>1) whichBmin<-whichBmin[1]
    if(y[whichBmin]==projYrs) return((B[whichBmin]-Btarget)^2)
    # project weekly over two years and take minimum
    res<-ode(y=res[whichBmin,c(2:ncol(res))],times=seq(0,2,0.02),func=fnFoodWeb_dBdt,parms=parms,method="rk4")
    return((min(res[,2])[1]-Btarget)^2)
    
  } else { # gamma too high and population to zero before projYrs
    return(  (B[nYrsGTE0]-(B[nYrsGTE0]-B[(nYrsGTE0+1)])*projYrs-Btarget)^2 )
  } # end else
}



