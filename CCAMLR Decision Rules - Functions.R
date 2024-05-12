# CCAMLR Decision Rules - Functions

# 1.0 General ####


# 2.0 Krill ####

#     2.1 Pella-Tomlinson 1969 ####

dP<-function(B,par){B*par$r/par$phi*(1-(B/par$K)^par$phi)}
  
fnPellaTomlinson_dBdt<-function(t,y,parms){
  catch<-if(length(parms$catch)>1) parms$catch[floor(t)+1] else parms$catch 
  list(dP(y,parms$PT$par)-fnM(y,parms$PT)-catch)
}


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



# 6.0 Plotting

fnPlotPellaTomlinson<-function(PT,PlotColours,HollingExamples=NULL,Krange=c(0,1),Ylim=NULL,plotMrate=TRUE,plotMtypeH=TRUE,catch=0){
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
    LineName<-"Rate"
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
      Line    <- rep(paste0(names(H)[h]),length(B))
      return(data.frame(Line,B,Change))
    },HollingExamples,B,PT)))           
  } # end if
  pd$Line<-factor(pd$Line,levels=unique(pd$Line))
  lineColours<-sapply(as.vector(levels(pd$Line)),function(l,PlotColours){if(l %in% names(PlotColours)) return(as.vector(PlotColours[l])) else 
    pc<-"black"; return(pc)},PlotColours)
  
  Maxima<-pd$B[pd$Line=="Productivity"][which(pd$Change[pd$Line=="Productivity"]==(max(pd$Change[pd$Line=="Productivity"])))]
  names(Maxima)="Productivity"
  usePD<-pd[!(pd$Line=="Productivity"),]
  usePD$Line<-factor(usePD$Line, levels=unique(usePD$Line))
  P0<-pd$Change[pd$Line=="Productivity"]
  Maxima<-c(Maxima,tapply(usePD,usePD$Line,function(pd,P0){pd$B[which((P0-pd$Change)==max(P0-pd$Change))]},P0))
  # names(Maxima)<-levels(pd$Line)
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
  p<-p+scale_colour_manual(values=as.vector(lineColours))
  p<- p+theme(legend.key=element_blank()) # removes boxes from around lines in legend
  
  # now do colours and line types as chosen                   
  #  ,colour="darkgreen",linetype="dashed")+geom_vline(xintercept=BatPmaxWithMfn,colour="darkgreen",linetype="dotted")
  p<-p+geom_vline(xintercept=Maxima,colour=as.vector(lineColours),linetype=rep("dotted",length(lineColours)))
  
  # annotations
  p<- p+annotate("text",x=PT$B0*1.05, y=(Ylim[2]-Ylim[1])*1E-1, label="B[0]",parse=TRUE,hjust=0,vjust=0)
  
  p<- p+annotate("text",x=Maxima[1], y=(Ylim[2]-Ylim[1])*5E-2, label="hat(P): M (0)"      ,parse=TRUE,hjust=0,vjust=0.5,angle=90,size=3)
  p<- p+annotate("text",x=Maxima[2], y=(Ylim[2]-Ylim[1])*5E-2, label="hat(P): M (rate)"   ,parse=TRUE,hjust=0,vjust=0.5,angle=90,size=3)
  p<- p+annotate("text",x=Maxima[3], y=(Ylim[2]-Ylim[1])*5E-2, label="hat(P): M (whale)"  ,parse=TRUE,hjust=0,vjust=0.5,angle=90,size=3)
  p<- p+annotate("text",x=Maxima[4], y=(Ylim[2]-Ylim[1])*5E-2, label="hat(P): M (penguin)",parse=TRUE,hjust=0,vjust=0.5,angle=90,size=3)
  
  return(p)
} # end fn


fnPlotPopWithGamma<-function(Gamma # vector with names for mortality type
                             ,parms,ProjYears,Xlim=NULL,XticksN=7,Ylim=NULL,YticksN=6,useStatus=TRUE, doLegend=TRUE,LegendTopLeft=c(0.6,0.25),doAxisLabels=TRUE,LineColours=NULL){
  
  LineColours<-plotColours

  # create dataframe ####
  y0<-parms$PT$B0
  names(y0)<-"B"
  plotYrs<-ProjYears$YearsTotal
  if(is.null(parms$catch)) parms$catch<-rep(0, ProjYears$YearsTotal)
  
  pd<-do.call(rbind,lapply(c(1:length(Gamma)),function(g,Gamma,parms,y0,plotYrs,yrsFishery){
    G<-Gamma[g]
    MortalityType<-names(G)
    H<-parms$HollingExamples
    yrsFishery<-yrsFishery[yrsFishery<=plotYrs] # ensure NAs do not enter the data
    if(MortalityType=="Rate"){
      parms$PT$par$useMrate<-TRUE
      parms$catch[yrsFishery]<-G*parms$PT$B0
      res<-as.data.frame(ode(y=y0,times=seq(0,plotYrs,0.2),func=fnPellaTomlinson_dBdt,parms=parms,method="rk4"))
      Mortality<-rep(MortalityType,nrow(res))
      Gamma<-rep(G,nrow(res))
      return(cbind(Mortality,Gamma,res))
    } else {
      parms$PT$par$Holling<-H[[MortalityType]]
      parms$PT$par$useMrate<-FALSE
      PB0<-parms$PT$B0*parms$PT$par$M
      parms$PT$par$Mmax<-nlm(findMmax,PB0,parms$PT,parms$PT$B0,PB0)$estimate
      parms$catch[yrsFishery]<-G*parms$PT$B0
      res<-as.data.frame(ode(y=y0,times=seq(0,plotYrs,0.2),func=fnPellaTomlinson_dBdt,parms=parms,method="rk4"))
      Mortality<-rep(MortalityType,nrow(res))
      Gamma<-rep(G,nrow(res))
      return(cbind(Mortality,Gamma,res))
    } # end else
  },Gamma,parms,y0,plotYrs,ProjYears$Fishery))
  
  pd$Mortality<-factor(pd$Mortality,levels=unique(pd$Mortality))
  
  # Plotting routine ####  
  if(useStatus)                pd$B<-pd$B/pd$B[1]
  if(is.null(Xlim))            Xlim<-c(0,max(pd$time,na.rm=TRUE))
  if(is.null(Ylim))            Ylim<-c(0,max(pd$B,na.rm=TRUE))
  
  pd<-pd[pd$time>=Xlim[1] & pd$time<=Xlim[2] & !is.na(pd$B),] # restrict plot to x limits & remove NAs from dataframe
  
  p<-ggplot(pd)+labs(x="Year",y="Status")
  p<-p+theme(panel.background = element_rect(fill="white",colour = "black")
             ,panel.grid.major = element_blank()
             ,panel.grid.minor = element_blank()
             ,axis.line=element_line(colour="black")
             ,axis.text = element_text(size=14)
             ,axis.title=element_text(size=18)
             ,legend.position = LegendTopLeft
             ,legend.background = element_rect(fill="white",size=0.5,linetype="solid",colour="black")
  )# end theme
  p<- p+scale_x_continuous(breaks=seq(Xlim[1],Xlim[2],length.out=XticksN)
                           ,limits=Xlim
                           ,expand=expansion(mult = c(0,0.05)))
  p<- p+scale_y_continuous(limits=Ylim
                           ,breaks=seq(Ylim[1],Ylim[2],length.out=YticksN)
                           ,expand=expansion(mult = c(0,0.05)))
  p<- p+guides(
    x = guide_axis(minor.ticks = TRUE)
    ,y = guide_axis(minor.ticks = TRUE)
  )
  
  p<- p+geom_hline(yintercept=c(0.2,0.75,1.0),colour=c("black","black","black"),linetype=c("dotted","dotted","dotted"))
  p<- p+geom_vline(xintercept=c((ProjYears$Fishery[1]-1)
                                ,ProjYears$Fishery[length(ProjYears$Fishery)]
                                ,(ProjYears$Fishery[length(ProjYears$Fishery)]+ProjYears$Recovery))
                   ,colour=c("black","black","black"),linetype=c("dashed","dashed","dashed"))
  
  p<- p+geom_line(aes(x=time,y=B,colour=Mortality),linewidth=1)
  if(!is.null(LineColours)) {
    lColours<-c("black",LineColours[c("Whale","Penguin")])
    names(lColours)<-NULL
    p<- p+scale_colour_manual(values=lColours, labels= c(bquote(paste( "Rate:      ",gamma," = ",.(round(Gamma[1],3))))
                                                         ,bquote(paste("Whale:    ",gamma," = ",.(round(Gamma[2],3))))
                                                         ,bquote(paste("Penguin: ",gamma," = ",.(round(Gamma[3],3))))
    ))
    #levels(pd$Mortality),' : ',expression(gamma),' = ',round(Gamma,3)))
  }
  p<- p+theme(legend.key=element_blank())
  return(p)
} # end fn

fnPlotFWwithGammaSinglePreds<-function(GammaFW # vector with names for mortality type
                             ,parms,ProjYears,ChangingK=FALSE
                             ,Xlim=NULL,XticksN=7,Ylim=NULL,YticksN=6,useStatus=TRUE
                             , doLegend=TRUE,LegendTopLeft=c(0.6,0.25),doAxisLabels=TRUE
                             ,LineColours=NULL){
  
  # create dataframe ####
 
  pd<-do.call(rbind,lapply(c(1:length(parms$pC)),function(p,parms,plotYrs,Gamma,yrsFishery,ChangingK,useStatus){
    Scenario<-names(parms$pC[p])
    par<-parms
    par$pC<-parms$pC[p]
    yrsFishery<-yrsFishery[yrsFishery<=plotYrs] # ensure NAs do not enter the data
    G<-Gamma[p]
    fwB0<-c(par$pP$B0,as.vector(sapply(par$pC,function(p){p$B0})))  # biomasses of krill/preds time 0 (from Ecopath)
    names(fwB0)<-c(par$pP$name,as.vector(sapply(par$pC,function(p){p$name})))

    # do projection

    par$catch<-rep(0,plotYrs)

    par$catch[yrsFishery]<-G*par$pP$B0  # set the catch for the period of the fishery 

    par$pP$modelK$changingK <-ChangingK
    
    res<-ode(y=fwB0,times=seq(0,plotYrs-1,0.2),func=fnFoodWeb_dBdt,parms=par,method="rk4")
    # develop dataframe with relevant grouping variables
   if(useStatus){
     # divide each biomass series by first in series
     r1<-do.call(rbind,lapply(seq(2,ncol(res),1),function(s,B,B0){
       return(data.frame(Scenario=rep(paste0(Scenario,"_",dimnames(B)[[2]][s]),nrow(B))
                         ,Gamma=rep(round(G,3),nrow(B))
                         ,Taxon=rep(dimnames(B)[[2]][s],nrow(B))
                         ,Year=B[,1]
                         ,Status=B[,s]/B0[(s-1)]))
     },res,fwB0))
   } else {
     r1<-do.call(rbind,lapply(seq(2,ncol(res),1),function(s,B,B0){
       return(data.frame(Scenario=rep(paste0(Scenario,"_",dimnames(B)[[2]][s]),nrow(B))
                         ,Gamma=rep(round(G,3),nrow(B))
                         ,Taxon=rep(dimnames(B)[[2]][s],nrow(B))
                         ,Year=B[,1]
                         ,Status=B[,s]))
     },res,fwB0))
   }
    return(r1)
    },parms,ProjYears$YearsTotal,GammaFW,ProjYears$Fishery,ChangingK,useStatus))

  pd$Scenario<-factor(pd$Scenario,levels=unique(pd$Scenario))

  # Plotting routine ####  
  if(is.null(Xlim))            Xlim<-c(0,max(pd$Year,na.rm=TRUE))
  if(is.null(Ylim))            Ylim<-c(0,max(pd$Status,na.rm=TRUE))
  
  pd<-pd[pd$Year>=Xlim[1] & pd$Year<=Xlim[2] & !is.na(pd$Status),] # restrict plot to x limits & remove NAs from dataframe
  pColTaxa<-as.vector(sapply(levels(pd$Scenario),function(s,pd){unique(pd$Taxon[pd$Scenario==s])},pd))

  p<-ggplot(pd)
  p<-p+theme(panel.background = element_rect(fill="white",colour = "black")
             ,panel.grid.major = element_blank()
             ,panel.grid.minor = element_blank()
             ,axis.line=element_line(colour="black")
             ,axis.text = element_text(size=14)
             ,axis.title=element_text(size=18)
             ,legend.position = LegendTopLeft
             ,legend.background = element_rect(fill="white",size=0.5,linetype="solid",colour="black")
  )# end theme
  p<- p+scale_x_continuous(breaks=seq(Xlim[1],Xlim[2],length.out=XticksN)
                           ,limits=Xlim
                           ,expand=expansion(mult = c(0,0.05)))
  p<- p+scale_y_continuous(limits=Ylim
                           ,breaks=seq(Ylim[1],Ylim[2],length.out=YticksN)
                           ,expand=expansion(mult = c(0,0.05)))
  p<- p+guides(
    x = guide_axis(minor.ticks = TRUE)
    ,y = guide_axis(minor.ticks = TRUE)
  )
  
  p<- p+geom_hline(yintercept=c(0.2,0.75,1.0),colour=c("black","black","black"),linetype=c("dotted","dotted","dotted"))
  p<- p+geom_vline(xintercept=c((ProjYears$Fishery[1]-1)
                                ,ProjYears$Fishery[length(ProjYears$Fishery)]
                                ,(ProjYears$Fishery[length(ProjYears$Fishery)]+ProjYears$Recovery))
                   ,colour=c("black","black","black"),linetype=c("dashed","dashed","dashed"))
  
  p<- p+geom_line(aes(x=Year,y=Status,col=Scenario,linetype=Scenario),linewidth=1)
  if(!is.null(LineColours)) {
    lColours<-LineColours[pColTaxa]
    names(lColours)<-NULL
    p<- p+scale_linetype_manual(values=c("solid","solid","twodash","twodash")
                                , labels= c("Krill (Whale scenario)","Whale","Krill (Penguin scenario)","Penguin"))
    
    p<- p+scale_colour_manual(values=lColours
                              , labels= c("Krill (Whale scenario)","Whale","Krill (Penguin scenario)","Penguin"))
    p<-p+labs(x="Year",y="Status"
              ,colour=bquote(paste( "Scenario (",gamma," = ",.(round(Gamma[1],3)),")"))
              ,linetype=bquote(paste( "Scenario (",gamma," = ",.(round(Gamma[1],3)),")")))  
    }
  p<- p+theme(legend.key=element_blank())
  return(p)
} # end fn



fnPlotFWwithGammaTwoPreds<-function(GammaFW # vector with names for mortality type
                                         ,parms,ProjYears,ChangingK=FALSE
                                         ,startWhalesX=1
                                         ,Xlim=NULL,XticksN=7,Ylim=NULL,YticksN=6,useStatus=TRUE
                                         , doLegend=TRUE,LegendTopLeft=c(0.6,0.25),doAxisLabels=TRUE
                                         ,LineColours=NULL){
   
    # create dataframe ####
    
      plotYrs<-ProjYears$YearsTotal
      yrsFishery<-ProjYears$Fishery
      yrsFishery<-yrsFishery[yrsFishery<=plotYrs] # ensure NAs do not enter the data
      G<-GammaFW
      fwB0<-c(parms$pP$B0,as.vector(sapply(parms$pC,function(p){p$B0})))  # biomasses of krill/preds time 0 (from Ecopath)
      names(fwB0)<-c(parms$pP$name,as.vector(sapply(parms$pC,function(p){p$name})))

      fwStart<-fwB0
      fwStart["Whale"]<-fwStart["Whale"]*startWhalesX

      # do projection
      
      parms$catch<-rep(0,plotYrs)
      parms$catch[yrsFishery]<-G*parms$pP$B0  # set the catch for the period of the fishery 
      
      parms$pP$modelK$changingK <-ChangingK
      
      res<-ode(y=fwStart,times=seq(0,plotYrs-1,0.2),func=fnFoodWeb_dBdt,parms=parms,method="rk4")
      # develop dataframe with relevant grouping variables
      if(useStatus){
        # divide each biomass series by first in series
        pd<-do.call(rbind,lapply(seq(2,ncol(res),1),function(s,B,B0,G){
          return(data.frame(Gamma=rep(round(G,3),nrow(B))
                            ,Taxon=rep(dimnames(B)[[2]][s],nrow(B))
                            ,Year=B[,1]
                            ,Status=B[,s]/B0[(s-1)]))
        },res,fwB0,G))
      } else {
        pd<-do.call(rbind,lapply(seq(2,ncol(res),1),function(s,B,G){
          return(data.frame(Gamma=rep(round(G,3),nrow(B))
                            ,Taxon=rep(dimnames(B)[[2]][s],nrow(B))
                            ,Year=B[,1]
                            ,Status=B[,s]))
        },res,fwB0,G))
      }

    pd$Taxon<-factor(pd$Taxon,levels=unique(pd$Taxon))
    
    # Plotting routine ####  
    if(is.null(Xlim))            Xlim<-c(0,max(pd$Year,na.rm=TRUE))
    if(is.null(Ylim))            Ylim<-c(0,max(pd$Status,na.rm=TRUE))
    
    pd<-pd[pd$Year>=Xlim[1] & pd$Year<=Xlim[2] & !is.na(pd$Status),] # restrict plot to x limits & remove NAs from dataframe
    pColTaxa<-names(fwB0)
    
    p<-ggplot(pd)
    p<-p+labs(x="Year",y="Status"
              ,colour=bquote(paste( "Krill catch ",gamma," = ",.(round(G,3)))))
    
    p<-p+theme(panel.background = element_rect(fill="white",colour = "black")
               ,panel.grid.major = element_blank()
               ,panel.grid.minor = element_blank()
               ,axis.line=element_line(colour="black")
               ,axis.text = element_text(size=14)
               ,axis.title=element_text(size=18)
               ,legend.position = LegendTopLeft
               ,legend.background = element_rect(fill="white",size=0.5,linetype="solid",colour="black")
    )# end theme
    p<- p+scale_x_continuous(breaks=seq(Xlim[1],Xlim[2],length.out=XticksN)
                             ,limits=Xlim
                             ,expand=expansion(mult = c(0,0.05)))
    p<- p+scale_y_continuous(limits=Ylim
                             ,breaks=seq(Ylim[1],Ylim[2],length.out=YticksN)
                             ,expand=expansion(mult = c(0,0.05)))
    p<- p+guides(
      x = guide_axis(minor.ticks = TRUE)
      ,y = guide_axis(minor.ticks = TRUE)
    )
    
    p<- p+geom_hline(yintercept=c(0.2,0.75,1.0),colour=c("black","black","black"),linetype=c("dotted","dotted","dotted"))
    p<- p+geom_vline(xintercept=c((ProjYears$Fishery[1]-1)
                                  ,ProjYears$Fishery[length(ProjYears$Fishery)]
                                  ,(ProjYears$Fishery[length(ProjYears$Fishery)]+ProjYears$Recovery))
                     ,colour=c("black","black","black"),linetype=c("dashed","dashed","dashed"))
    
    p<- p+geom_line(aes(x=Year,y=Status,col=Taxon),linewidth=1)
    if(!is.null(LineColours)) {
      lColours<-LineColours[pColTaxa]
      names(lColours)<-NULL
      p<- p+scale_colour_manual(values=lColours
                                , labels= c("Krill","Whale","Penguin"))
         } # end if
    p<- p+theme(legend.key=element_blank())
    return(p)
  } # end fn
