fPathIn<-"./IO in/"
fPathOut<-"./IO out/"
load(paste0(fPathIn,"Fe_means.Rdata")) # df Iron concentration rows(MEASO area) cols(n, name, value)
load(paste0(fPathIn,"Si_means.Rdata")) # df Silicate concentration rows(MEASO area) cols(n, name, value)
load(paste0(fPathIn,"Insol_means.Rdata"))  # matrix Insolation W.m-2.d-1 rows(months) cols(MEASO area)
load(paste0(fPathIn,"MLD_means.Rdata"))  # matrix Mixed Layer Depth m rows(months) cols(MEASO area)
load(paste0(fPathIn,"CICE_means.Rdata"))  # matrix Sea ice percent cover % rows(months) cols(MEASO area) 
load(paste0(fPathIn,"Temp_MLD_means.Rdata"))  # matrix Temperature mixed layer oC rows(months) cols(MEASO area) 
load(paste0(fPathIn,"Temp_Deep_means.Rdata"))  # matrix Temperature MLD to 1000m oC rows(months) cols(MEASO area) 

TargetArea<-"AOA"
nFe_0<-Fe_means[Fe_means[,"MEASO_area"]==TargetArea,3]
nSi_0<-Si_means[Si_means[,"MEASO_area"]==TargetArea,3]
eInsol<-Insol_means[,TargetArea]
eMLD<-MLD_means[,TargetArea]
eCICE<-CICE_means[,TargetArea]
eTemp_MLD<-Temp_MLD_means[,TargetArea]
eTemp_Deep<-Temp_Deep_means[,TargetArea]



## RHS of NPZD system (from Melbourne-Thomas et al 2015)

fnSOwebDE <- function(k # vector element to read  (not used by JMT)
              ,x # X vector - NPZD
              ,u # mortality rate
              ,a # parameters
              ,mld # mixed layer depth for time step
              ,no3e # nitrate concentration for time step
              ,w # change in mixed layer depth for time step
              ,J # phytoplankton growth rate
  ## Grazing
  G <- (a[9]*a[10]*x[2]^2)/(a[9] + a[10]*x[2]^2)
  wp <- max(w,0)
  c(## Nitrate equation
    a[13]*x[4]+a[12]*x[3]-J*x[2]+(wp*(no3e-x[1]))/mld,
    ## Phytoplankton equation
    (J-wp/mld)*x[2]-a[6]*x[2]-a[7]*x[2]^2-G*x[3],
    ## Zooplankton equation
    (a[8]*G-a[12]-w/mld)*x[3]-a[11]*x[3]^2,
    ## Detritus equation
    (1-a[8])*G*x[3]+a[11]*x[3]^2+a[6]*x[2]+a[7]*x[2]^2-(a[13]+(wp+a[14])/mld)*x[4])
}

