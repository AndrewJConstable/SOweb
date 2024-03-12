rootEst<-"/Users/acon/Desktop/_w/_r/SOweb/EcopathEst_Old/"
saveEst<-"/Users/acon/Desktop/_w/_r/SOweb/EcopathEst/"
rootDiet<-"/Users/acon/Desktop/_w/_p/SOfoodweb & EcoVelocity/Ecosystem network synthesis/Ecopath files/Rinputs/"
fNames<-c("Ballerini2014","Ballerini2014_Stand","Dahood2019"
         ,"Gurney2014_Yr2000","Gurney2014_Stand","Hill2021"
         ,"Hill2021_Stand","Maldonado2016_Yr1900","Maldonado2016_Yr2008"
         ,"McCormack2019","Pinkerton2010","Pinkerton2010_Stand")

units<-365

sapply(fNames,function(f){
  print(f)
eMod<-fnReadEcopath(rootDiet,f)
QB<-eMod$params[,"Q.B"]
Diet<-eMod$diet
fEcopathEst<-paste0("EcopathEst_",f,"_day_OLS.Rdata")
load(paste0(rootEst,fEcopathEst))
nEst<-EcopathEst$PoolN
vEst<-EcopathEst$v_matrix
if(nrow(vEst)==nEst) saveRDS(EcopathEst,paste0(saveEst,fEcopathEst)) else {
 # update v matrix
    vNew<-matrix(0,ncol=nEst,nrow=nEst)
   # find columns where only one prey and delete first 8 elements
    r1<-which(apply(Diet,2,function(d){if(sum(d>0)==1) return(TRUE) else return(FALSE)}))
    
    # update vulnerability to 1 (because the solution was any vulnerability) and update ingesiton rate
    # from Q.B
    
    vNew[,r1]<-do.call(cbind,lapply(r1,function(j,v,n){
      v<-v[c(9:(nEst+8)),j]
      v[which(v>0)]<-1
      return(v)
      },vEst,nEst))
    EcopathEst$IngestionMax[r1]<-unlist(sapply(r1,function(j){
        print(list(which(vNew[,j]>0)[1],eMod$params[,"B"],vNew[,j],QB[j]/units))
             res<-fnCalcI(which(vNew[,j]>0)[1],eMod$params[,"B"],vNew[,j],QB[j]/units)
            # print(res)
             if(res>0) return(res) else return(QB[j]/units)}))
        # for other columns, delete last 8 elements
   r2<-c(1:nEst)[-r1]
   vNew[,r2]<-vEst[c(1:nEst),r2]
    EcopathEst$v_matrix<-vNew
    saveRDS(EcopathEst,paste0(saveEst,fEcopathEst))
    }  # end else
})
