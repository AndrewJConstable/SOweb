rootEcopathEst<-"/Users/acon/Desktop/_w/_r/SOweb/EcopathEst/20240304 original est/"
outEcopathEst<-"/Users/acon/Desktop/_w/_r/SOweb/EcopathEst/"
rootEcopath<-"/Users/acon/Desktop/_w/_p/SOfoodweb & EcoVelocity/Ecosystem network synthesis/Ecopath files/Rinputs/"
rootIOin<-"/Users/acon/Desktop/_w/_r/SOweb/IO in/"

fSOwebParams<-"SOwebPoolParams.csv"
SOwebParams<-read.csv(paste0(rootIOin,fSOwebParams))

SOweb_pools<-SOwebParams$VarName


Ecopath_Mods<-c("Dahood2019","Gurney2014_Yr2000","Hill2021"
                ,"Maldonado2016_Yr1900","Maldonado2016_Yr2008"
                ,"McCormack2019","Subramanium2020","Surma2014")
# note - the following do not need adjusting - see Hill et al 2021
# "Ballerini2014","Ballerini2014_Stand","Gurney2014_Stand","Hill2021_Stand"
# ,"Pinkerton2010","Pinkerton2010_Stand"


sapply(Ecopath_Mods,function(f){  # start sapply ecopath models
print(f)
eMod<-fnReadEcopath(rootEcopath,f)

B<-eMod$params[,"B"]
QBw<-eMod$params[,"Q.B"]
Diet<-eMod$diet
CF<-unlist(lapply(eMod$names[,"ConvSOvar"],function(p,SOw){
  mask<-SOw[,"VarName"] %in% p
  if(sum(mask)>0) return(SOw[which(mask),"gCfromWW"]) else return(0)
},SOwebParams))

QBc<-unlist(sapply(c(1:nrow(eMod$names)),function(j,B,QBw,d,CF){
  return(QBw[j]*sum(d[,j]*CF)/CF[j])
},B,QBw,Diet,CF))

fEcopathEst<-paste0("EcopathEst_",f,"_day_OLS.Rdata")
load(paste0(rootEcopathEst,fEcopathEst))

EcopathEst$IngestionMax<-EcopathEst$IngestionMax*QBc/QBw
EcopathEst$IngestionMax[is.na(EcopathEst$IngestionMax)]<-0

save(EcopathEst,file=paste0(outEcopathEst,fEcopathEst))
}) # end sapply ecopath models
