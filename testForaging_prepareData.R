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

Weddell Sea quantitative network model
https://github.com/EcoComplex/WeddellSea

###########################
# Southern Ocean Dietary database

# downloaded 19 February 2024
library(sohungry)
library(solong)
library(dplyr)

root<-"/Users/acon/Desktop/_w/_d/SCAR_Diet_Energetics/"
fSO_diet<-"scar_diet.csv"
dfSO_diet<-read.csv(paste0(root,fSO_diet))

taxaFamily<-unique(dfSO_diet[,c("predator_worms_class","predator_worms_order","predator_worms_family")])
taxaFamily<-taxaFamily[order(taxaFamily[,"predator_worms_class"]
                             ,taxaFamily[,"predator_worms_order"]
                             ,taxaFamily[,"predator_worms_family"])
                             ,]



# Southern Ocean dietary database
https://zenodo.org/records/7796465
Citation:
Scientific Committee on Antarctic Research. (2023). SCAR Southern Ocean Diet and Energetics Database (2023-04-04) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.7796465

variable name details:
  
  
# R libraries - sohungry and so

Github: https://github.com/SCAR/sohungry


options(repos = c(scar = "https://scar.r-universe.dev", CRAN = "https://cloud.r-project.org"))
install.packages("sohungry")

## or install from github
## install.packages("remotes") ## if needed
remotes::install_github("SCAR/sohungry")

install.packages("devtools")
library(devtools)
install_github("SCAR/solong")

taxaFamily<-unique(dfSO_diet[,c("predator_worms_class","predator_worms_order","predator_worms_family")])
taxaFamily<-taxaFamily[order(taxaFamily[,"predator_worms_class"]
                             ,taxaFamily[,"predator_worms_order"]
                             ,taxaFamily[,"predator_worms_family"])
                             ,]

# Predators - Birds

# Class = Aves, Order = Procellariformes (flying birds), Sphenisciformes (penguins)

# Flying birds - all Procellariiformes
# Small ocean penguins - Pygoscelis papua (gentoo), Eudyptes chrysolophus (macaroni) 
# Small ice penguins - Pygoscelis adeliae (Adelie)
# large ocean penguins - Aptenodytes patagonicus (King)
# large ice penguins - Aptenodytes forsteri

aFB<-list(level="predator_worms_genus",wormNames = c("Diomedea","Thalassarche","Phoebetria")) # family c("Diomedeidae"))
aPSO<-list(level="predator_worms_species",wormNames = c("papua","chrysolophus"))
aPSI<-list(level="predator_worms_species",wormNames = c("adeliae"))
aPLO<-list(level="predator_worms_species",wormNames = c("patagonicus"))
aPLI<-list(level="predator_worms_species",wormNames = c("forsteri"))

unique(tmp[,c("predator_name","predator_worms_order"
              ,"predator_worms_family","predator_worms_genus")])
taxaAves<-dfSO_diet[(dfSO_diet[,"predator_worms_order"]=="Procellariiformes" | 
                      dfSO_diet[,"predator_worms_genus"]=="Pygoscelis" |
                      dfSO_diet[,"predator_worms_genus"]=="Aptenodytes" |
                      dfSO_diet[,"predator_worms_genus"]=="Eudyptes"),c("predator_name","predator_worms_order"
                            ,"predator_worms_family","predator_worms_genus")]

taxaAves<-unique(taxaAves[c("predator_name","predator_worms_order"
                            ,"predator_worms_family","predator_worms_genus")])

black-browed albatross (Thalassarche melanophris)
grey-headed albatross (Thalassarche chrysostoma),
light-mantled sooty albatross (Phoebetria palpebrata)
wandering albatross (Diomedea exulans)

# simplify - do family Diomedeidae (covers the main species of albatross)

aFBd<-dfSO_diet[dfSO_diet[,aFB[["level"]]]%in%aFB[["wormNames"]],]
unique(aFBd[,c("prey_worms_class")]) #"prey_size_mean","prey_size_units"

unique(aFBd[aFBd[,"prey_worms_class"]=="Teleostei",c("predator_worms_genus","prey_name","prey_worms_family")]) #"prey_size_mean","prey_size_units"

aves_prey_classes<-c("Cephalopoda","Petromyzonti","Actinopterygii"
,"Malacostraca","Hexanauplia","Aves","Bivalvia"      
,"Polychaeta","Scyphozoa","Teleostei","Elasmobranchii"
,"Gastropoda","Anthozoa","Hydrozoa")
UsePreyNames<-aves_prey_classes[!(aves_prey_classes%in%c("Cephalopoda","Teleostei"))]
hist(aFBd[aFBd[,"prey_worms_class"]%in%UsePreyNames
            ,"prey_size_mean"])
unique(aFBd[aFBd[,"prey_worms_class"]%in%UsePreyNames & aFBd[,"prey_size_mean"]>100
          ,"prey_worms_class"])
# note that these are all benthic species (including less than 100. Therefor excluded as inshore or from fishery


UsePreyNames<-"Cephalopoda"
hist(aFBd[aFBd[,"prey_worms_class"]%in%UsePreyNames
          ,"fraction_occurrence"])

CephalopodFamilies <- unique(aFBd[aFBd[,"prey_worms_class"]%in%UsePreyNames
          ,"prey_worms_family"])

res<-lapply(CephalopodFamilies,function(f,d){
     df<-d[d[,"prey_worms_family"]==f,]
     N<-nrow(df)
     Size_min<-mean(df[,"prey_size_mean"],na.rm=TRUE)
     Size_max<-mean(df[,"prey_size_mean"],na.rm=TRUE)
     Diet_min<-min(df[,"fraction_occurrence"],na.rm=TRUE)
     Diet_median<-median(df[,"fraction_occurrence"],na.rm=TRUE)
     Diet_max<-max(df[,"fraction_occurrence"],na.rm=TRUE)
     return(list(N = N
                ,Size_min = Size_min
                ,Size_max = Size_max
                ,Diet_min = Diet_min
                ,Diet_median = Diet_median
                ,Diet_max = Diet_max))
},aFBd[aFBd[,"prey_worms_class"]%in%UsePreyNames,])
res<-data.frame(matrix(as.numeric(do.call(rbind,res)),ncol=6))
res<-cbind(CephalopodFamilies,res)
names(res)<-c("Family","N","Size_min","Size_max","Diet_min","Diet_median","Diet_max")
res<-res[order(res[,"N"],res[,"Diet_median"]),]



