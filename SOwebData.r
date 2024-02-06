# SOwebData.r
# Input data and parameters
#    1.1. Minimisation inputs ####
# time0 relates to 1st element of a temporal vector.  The matrix of state variables 
nYrs<-1
nStepsPerYr<-365 
Date1<-as.Date("2015-07-01")


#    1.2. Input Parameters ####

a<-list(  Nu = NULL
         ,Ph = list(Di = list( X     = 1 # position in state variable vector
                              ,MuMax = 1.44 # from Jeffery
                              ,alpha = 0.16    # (W m-2 d)-1 Hauck - note alpha*q{Chl:N}*q{N:C}
                              ,k_Fe  = 0.12    # micromole Fe m-3 Hauck
                              ,k_Si  = 4.0    # millimole Si m-3 Hauck
                              ,r_FeC = 0.005 # micromole Fe (millimole C)-1 Hauck
                              ,r_SiC = 0.8   # mole Si (mole C)-1 Hauck
         ) # 
                   ,Sm = list( X     = 2
                              ,MuMax = 0.66    # Jeffery
                              ,alpha = 0.088   # Hauck
                              ,k_Fe  = 0.02    # micromole Fe m-3 Hauck
                              ,r_FeC = 0.005   # micromole Fe (millimole C)-1 Hauck
                               )) #
         ,Zo = NULL
         ,De = NULL
       )

#    1.3. Constants

cE<-list(Lat = eLat
         ,par = list(w = 0.43, si = 0.1) # proportion of insolation available as PAR without (w) and with sea ice present (si) sea ice has attenuation of 1.5 m-1
         ,kw  = 0.04                     # light attenuation (note that this is converted to 1/0.04 for computational efficiency)
         ,Fe  = nFe_0                    # iron concentration in deep water
         ,Si  = nSi_0                    # silicic acid concentration in deep water
) # end eC


#    1.4. Input Environmental Data ####

fPathIn<-"./IO in/"
fPathOut<-"./IO out/"
load(paste0(fPathIn,"Fe_means.Rdata")) # df Iron concentration rows(MEASO area) cols(n, name, value)
load(paste0(fPathIn,"Si_means.Rdata")) # df Silicate concentration rows(MEASO area) cols(n, name, value)
load(paste0(fPathIn,"Lat_means.Rdata")) # df Latitude mean for MEASO areas rows(MEASO area) cols(n, name, value)
load(paste0(fPathIn,"Insol_means.Rdata"))  # matrix Insolation W.m-2.d-1 rows(months) cols(MEASO area)
load(paste0(fPathIn,"MLD_means.Rdata"))  # matrix Mixed Layer Depth m rows(months) cols(MEASO area)
load(paste0(fPathIn,"CICE_means.Rdata"))  # matrix Sea ice percent cover % rows(months) cols(MEASO area) 
load(paste0(fPathIn,"Temp_MLD_means.Rdata"))  # matrix Temperature mixed layer oC rows(months) cols(MEASO area) 
load(paste0(fPathIn,"Temp_Deep_means.Rdata"))  # matrix Temperature MLD to 1000m oC rows(months) cols(MEASO area) 


# 2.1 Managing time ####
DaysInMonth<-c(31,28,31,30,31,30,31,31,30,31,30,31)
DayMidMonth<-c(15,14,15,15,15,15,15,15,15,15,15,15)

