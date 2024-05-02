# Script for processing Ecopath files

library(ggplot2)
# Data inputs

root<-"/Users/acon/Desktop/_w/_p/SOfoodweb & EcoVelocity/Ecosystem network synthesis/Ecopath files/Rinputs/"

Ecopath_Mods<-list(WPA   = "Pinkerton2010"
                  ,EPA.APS = "Ballerini2014"
                  ,EPA.APN = "Dahood2019"
                  ,AOS   = "Hill2021"
                  ,CIA   = "McCormack2019"
                  ,CIS   = "Subramanium2020"
                  ,CIN   = "Gurney2014_Yr2000"
) # end Ecopath_Mods
                   

Ecopath_Mods_Stand<-list( WPA = "Pinkerton2010_Stand"
                         ,EPA = "Ballerini2014_Stand"
                         ,AOS = "Hill2021_Stand"
                         ,CIN = "Gurney2014_Stand"
)# end standardised models

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

# collapse all nodes into functional groups, give missing functional groups NA

fnEcopathFnGroupsBiomass<-function(ep,Grps,cg){
   r<-do.call(rbind,lapply(Grps,function(g,ep,cg){
      dg<-cg[cg[,"Group"] %in% g,"FnGroup"]
      r1<- c(g,sum(ep$params$B[ep$names$FnGroup %in% dg],na.rm=TRUE))
      names(r1)<-c("Group","Biomass")
      return(r1)
      },ep,cg))
   return(r)
    } #


fnEcopathFnGroupsDetritalPool<-function(ep,Grps,cg){
  NPP<-{
    useG<-ep$names$FnGroup %in% cg[cg[,"Group"] %in% c(1),"FnGroup"]
    sum(ep$params$B[useG]*ep$params$P.B[useG],na.rm=TRUE)
  } # end NPP
  res<-do.call(rbind,lapply(Grps,function(g,ep,cg){
               useG<-ep$names$FnGroup %in% cg[cg[,"Group"] %in% g,"FnGroup"]
               G<-ep$names$N[useG]
               r<-do.call(rbind,lapply(G,function(g,ep){
                 useg<-ep$names$N %in% g
                 UA<-if(!is.na(ep$params$Q.B[useg])) ep$params$B[useg]*ep$params$Q.B[useg]*ep$params$UA[useg] else NA
                 OM<-ep$params$B[useg]*ep$params$P.B[useg]*(1-ep$params$EE[useg])
                 R <- if(!is.na(ep$params$Q.B[useg]) & ep$params$Q.B[useg]>0) ep$params$B[useg]*(ep$params$Q.B[useg]*(1-ep$params$UA[useg])-ep$params$P.B[useg]) else NA
                 return(data.frame(UA=UA,OM=OM,R=R))
               },ep))
               
               det<-ep$names$FnGroup %in% cg[cg[,"Group"] %in% c(8),"FnGroup"]
               nrowDet<-sum(det)
               propQ<-rep(0,sum(useG)); DC<-rep(NA,sum(useG))
               if(nrowDet>0){
                   nrowDet<-sum(det)
                   if(sum(useG)>0){
                     dietDet<-as.matrix(ep$diet[det,useG],nrow=nrowDet)
                     propQ<-apply(dietDet,2,sum,na.rm=TRUE)
                     DC<-ep$params$B[useG]*ep$params$Q.B[useG]*propQ
                   }}
          r1<-c(g,sum(r$UA,na.rm=TRUE)
                  ,sum(r$OM,na.rm=TRUE)
                  ,sum(DC,na.rm=TRUE)
                  ,sum(r$R,na.rm=TRUE)
                  ,sum(r$UA*(1-propQ),na.rm=TRUE)
                  ,sum(r$OM*(1-propQ),na.rm=TRUE)
                  ,sum(r$R*(1-propQ),na.rm=TRUE)
          )    
              names(r1)<-c("Group","UA","OM","DC","R","UApp","OMpp","Rpp")
    return(r1)
  },ep,cg))

  return(list(NPP = NPP, Data = res))
} 

# Primary production
# select all FnGroup==1, return total productivity as sum(P.B1 x B1)

fnEcopathPrimaryProduction<-function(ep){
PP<-ep$names$FnGroup==1
d<-ep$params[PP,c("B","P.B")]
return(sum(d[,"B"]*d[,"P.B"],na.rm=TRUE))
} # end function

# Summarise functional groups relative to total primary production
# B/P1  

##############################################################################
# 1. read standard parameters ####
eStand<-fnReadEcopath(root,Ecopath_Mods_Stand[4])
eStand$params


# Functional group definitions and combinations for plotting ####
# what FnGroups are represented

#            1 = primary producers
#            2 = zooplankton
#            3 = pelagic fish and squid
#            4 = birds and marine mammals (not including apex predators)
#            5 = apex predators
#            6 = benthic secondary producers
#            7 = benthic predators
#            8 = bentho-pelagic fish and squid
#            9 = detrital pool
#           10 = bacteria
CombinedGroups<-matrix(
  c(    1,1 # primary producers
      , 2,3 # zooplankton
      , 3,4 # pelagic fish and squid
      , 4,5 # birds and marine mammals
      , 5,5 # birds and marine mammals
      , 6,6 # benthos
      , 7,7 # bentho-pelagic fish and squid
      , 8,7 # bentho-pelagic fish and squid
      , 9,8 # detritus
      ,10,2 # bacteria
  ),ncol=2,byrow=TRUE) # end dataframe
dimnames(CombinedGroups)[[2]]<-c("FnGroup","Group")

GroupNames<-c("Protists","Bacteria","Zooplankton","Pelagic Fish & Squid","Birds & Marine Mammals","Benthos","Bentho-pelagic Fish & Squid","Detritus")
##############################################################################
# 2. process Ecopath models ####

fnBiomassRelPP<-function(epMods){
  mods<-names(epMods)
  do.call(rbind,lapply(mods,function(m,epMods){
    ep<-fnReadEcopath(root,epMods[[m]])
    Groups<-seq(1,length(GroupNames),1) # sort(unique(CombinedGroups[CombinedGroups[,"FnGroup"] %in% unique(ep$names$FnGroup),"Group"]))
    b<-as.data.frame(fnEcopathFnGroupsBiomass(ep,Groups,CombinedGroups))
    b[,"Biomass"]<-b[,"Biomass"]/b[b[,"Group"]==1,"Biomass"]
    b$Blog10<-log10(b[,"Biomass"])
    b$Blog10[is.infinite(b$Blog10)]<-NA
    Model<-rep(m,nrow(b))
    return(cbind(Model,b))
  },epMods))
} # end fnBiomassRelPP

fnDetritalPoolRelNPP<-function(epMods){
  mods<-names(epMods)
  r<-lapply(mods,function(m,epMods){
    ep<-fnReadEcopath(root,epMods[[m]])
    Groups<-seq(1,length(GroupNames),1) # sort(unique(CombinedGroups[CombinedGroups[,"FnGroup"] %in% unique(ep$names$FnGroup),"Group"]))
    r<-fnEcopathFnGroupsDetritalPool(ep,Groups,CombinedGroups) # returns list NPP (Net Primary Productivity), UA (Unassimilated), OM (Other Mortality) = (1-EE)
    dp<-as.data.frame(r$Data)
    Energy<-data.frame(NPP=r$NPP,Detritus=sum(r$Data[,"DC"],na.rm=TRUE))
    
    dp$RelDetPool<-(r$Data[,"UA"]+r$Data[,"OM"])/(r$NPP+sum(r$Data[,"DC"],na.rm=TRUE)) # hard to disentangle detritus as new energy
    dp$RDPlog10<-log10(dp[,"RelDetPool"])
    dp$RDPlog10[is.infinite(dp$RDPlog10)]<-NA
    Model<-rep(m,nrow(dp))
    dp<-cbind(Model,dp)
    Model<-rep(m,nrow(Energy))
    Energy<-cbind(Model,Energy)
    return(list(dp=dp,Energy=Energy))
  },epMods)
  return(list(dp     =do.call(rbind,lapply(r,function(m){m$dp}))
              ,Energy=do.call(rbind,lapply(r,function(m){m$Energy}))))
} # end fnDetritalPoolRelNPP


# 2. Biomass of trophic levels relative to primary producers ####



#Brel<-fnBiomassRelPP(Ecopath_Mods_Stand)
Brel<-fnBiomassRelPP(Ecopath_Mods)
Brel[,"Model"]<-factor(Brel[,"Model"],levels = unique(Brel[,"Model"]))
Brel$Taxa<-factor(GroupNames[Brel$Group],levels=GroupNames)
Brel$plot0<-Brel[,"Blog10"]
Brel$plot0[Brel$plot0>0.01 | Brel$plot0<(-0.01)]<-NA

pData<-Brel[Brel$Group<8,] # leave out detritus
yAdj<- (-3)
yAdjMin<-0
yAdjMax<-4.5
yMin<-(-3)
yMax<-yMin+yAdjMax-yAdjMin
yInt<-1
p<-ggplot(pData,aes(x=Group,y=Blog10-yAdj,fill=Taxa))+geom_bar(stat="identity", position=position_dodge()) 
p<-p + scale_fill_viridis_d()   #+ scale_fill_manual(values=c("darkgreen","brown1","orange","blue","grey","brown","darkblue","black"))
#p<-p+geom_point(aes(x=Group,y=plot0,colour=Taxa))  + scale_color_manual(values=c("darkgreen","brown1","orange","blue","grey","brown","darkblue","black"))
p<-p+labs(y="Biomass relative to Protists (Log10 )")
p<-p+theme(axis.title.x=element_blank()
           ,axis.text.x=element_blank()
           ,axis.ticks.x=element_blank()
           ) # end theme
p<-p+scale_y_continuous(limits=c(yAdjMin,yAdjMax),breaks=seq(yAdjMin,(yAdjMax),yInt), labels=seq(yMin,yMax,yInt))
p<-p+facet_grid(~Model)
p



# 3. Respiration, unassimilated food and (1-EE) relative to primary production ####

#    Respiration and Contributions to detrital pool of higher trophic levels
#        relative to primary produciton (B * P/B) 
#        Categories: sum[Respiration, Unassimilated food, (1-EE)]  = (Primary production + detritus used)
#        output detrital contribution to production 
#        detrital pool (carcass + detritus) = primary production less respiration
#        

DetritusPool<-fnDetritalPoolRelNPP(Ecopath_Mods)
DPrel<-DetritusPool$dp
DPrel[,"Model"]<-factor(DPrel[,"Model"],levels = unique(DPrel[,"Model"]))
DPrel$Taxa<-factor(GroupNames[DPrel$Group],levels=GroupNames)
DPrel$plot0<-DPrel[,"RDPlog10"]
DPrel$plot0[DPrel$plot0>0.01 | DPrel$plot0<(-0.01)]<-NA

pData<-DPrel[DPrel$Group<8,] # leave out detritus
yAdj<- (-5)
yAdjMin<-0
yAdjMax<-5
yMin<-(-5)
yMax<-yMin+yAdjMax-yAdjMin
yInt<-1
p<-ggplot(pData,aes(x=Group,y=RDPlog10-yAdj,fill=Taxa))
p<-p+geom_bar(stat="identity", position=position_dodge())
p<-p + scale_fill_viridis_d()   #+ scale_fill_manual(values=c("darkgreen","brown1","orange","blue","grey","brown","darkblue","black"))
p<-p+labs(y="Proportion Energy Input to Detrital Pools (Log10 )")
p<-p+theme(axis.title.x=element_blank()
           ,axis.text.x=element_blank()
           ,axis.ticks.x=element_blank()
) # end theme
p<-p+scale_y_continuous(limits=c(yAdjMin,yAdjMax),breaks=seq(yAdjMin,(yAdjMax),yInt), labels=seq(yMin,yMax,yInt))
p<-p+facet_grid(~Model)
p


#########################################
# nMDS of food web models
library(vegan)
library(ggrepel)
library(grid)

Model<-names(Ecopath_Mods)
d<-do.call(rbind,lapply(Model,function(m,d){return(d[d[,"Model"]==m,"Biomass"])},Brel[Brel$Group<8,c("Model","Biomass")]))
d
dnmds<-d[,-c(1,2)]
dnmds
set.seed(123)
rnmds=metaMDS(dnmds,distance="bray")
rnmds
snmds<-as.data.frame(scores(rnmds)$sites)
snmds<-cbind(Model,snmds)


p <- ggplot(snmds, aes(x = NMDS1, y = NMDS2, label = Model)) 
p <- p+  geom_point(size = 5, aes(colour = Model))
p <- p+ geom_text_repel()
p <- p+ geom_text(x=0.5, y=0.4, label="Stress < 0.001")
p <- p+ theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold")
        , axis.text.x = element_text(colour = "black", face = "bold", size = 12)
        ,legend.position="none"
#       ,legend.text = element_text(size = 12, face ="bold", colour ="black") 
#       ,legend.position = "right"
#        ,legend.title = element_text(size = 14, colour = "black", face = "bold")
        , axis.title.y = element_text(face = "bold", size = 14) 
        , axis.title.x = element_text(face = "bold", size = 14, colour = "black") 
        , panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2)
      #  ,legend.key=element_blank()
        )# end theme
p <- p+ labs(x = "NMDS1", y = "NMDS2")

p

ggsave("NMDS.svg")

