# SOwebData.r
# Input data and parameters
#    1.1. Minimisation inputs ####
# time0 relates to 1st element of a temporal vector.  The matrix of state variables 
nYrs<-1
nStepsPerYr<-365 
Date1<-as.Date("2015-07-01")


#    1.2. Input Parameters ####

a<-list( # list of taxa with their parameters and functions as needed (first letter of name is group, next 2 letters are taxon, 4th & more letter relates to subpool)
       
  nFeM =  list(Name    = "Nut-Iron-MLD"
              ,X0       = 20 # possible starting value
              ,Attr    = list( Units = "umole.m-3" # micromole
                              ,Which_C_ratio = "r_FeC"
                          )  
              ,Consume = # from detrital pool
                list(params = list(NULL
                                  ) # end list
                    ,fn     = "")
              ,Produce  = NULL # no production of nutrients
              ,Detritus = NULL # no function here
              ,Import   = NULL # from change in MLD and sea ice
              ,Export   = NULL # from change in MLD and sea ice and scavenging
  ) # end nFeM
  ,nFeSI =  list(Name    = "Nut-Iron-SeaIce"
                 ,X0       = 20 # possible starting value
                 ,Attr    = list(Units="umol.m-3"
                                 ,Which_C_ratio = "r_FeC"
                 )        
                 ,Consume = # for phytoplankton, consumption requires an estimate of production (the growth function below simply requires a conversion to carbon)
                   list(params = list(NULL
                   ) # end list
                   ,fn     = "")
                 ,Produce  = NULL # no production of nutrients
                 ,Detritus = NULL # no function here
                 ,Import   = NULL # from change in MLD and sea ice
                 ,Export   = NULL # from change in MLD and sea ice and scavenging
  ) # end nFeSI
  ,nSiM =  list(Name    = "Nut-Silicate-MLD"
               ,X0       = 20 # possible starting value
               ,Attr    = list( Units = "mmole.m-3" # millimole
                                ,Which_C_ratio = "r_SiC"
               )  
               ,Consume = # from detrital pool
                 list(params = list(NULL
                 ) # end list
                 ,fn     = "")
               ,Produce  = NULL # no production of nutrients
               ,Detritus = NULL # no function here
               ,Import   = NULL # from change in MLD and sea ice
               ,Export   = NULL # from change in MLD and sea ice and scavenging
  ) # end nSiM
 ,nSiSI =  list(Name    = "Nut-Silicate-Sea ice"
               ,X0       = 20 # possible starting value
               ,Attr    = list( Units = "mmole.m-3" # millimole
                                ,Which_C_ratio = "r_SiC"
               )  
               ,Consume = # from detrital pool
                 list(params = list(NULL
                 ) # end list
                 ,fn     = "")
               ,Produce  = NULL # no production of nutrients
               ,Detritus = NULL # no function here
               ,Import   = NULL # from change in MLD and sea ice
               ,Export   = NULL # from change in MLD and sea ice and scavenging
 ) # end nSiSI

   
  ,pDi = list( Name    = "Phyto-diatom" 
              ,X0       = 30 # possible starting value
              ,Attr   = list( Units = "g.m-2"
                             ,MassToMole = 1/12  # ensure correct units 1 g : 12 mole (if kg then = 1/12000 etc)
                             ,WtC   = NULL     # g.Carbon individual weight
                             ,WW_C  = NULL     # g.WetWeight per carbon (g.carbon)-1
                             ,r_FeC = 0.005   # micromole Fe (millimole C)-1 Hauck
                             ,r_SiC = 0.8     # mole Si (mole C)-1 Hauck
                             ) # end Attributes
              ,Consume = # for phytoplankton, consumption requires an estimate of production (the growth function below simply requires a conversion to carbon)
                              list(params = list(food   = c("nFeM","nSiM") # names of taxa being consumed - from names(a)
                                                ,MuMax = 1.44 # from Jeffery
                                                ,alpha = 0.16            # (W m-2 d)-1 Hauck - note alpha*q{Chl:N}*q{N:C}
                                                ,k    =  c(0.12,4.0)    # micromole.m-3, millimole.m-3 Hauck  - vector of half saturation constants for each food
                                                 ) # end list
                                 ,fn     = "fnPhConsume")
              ,Produce  = NULL # no production of nutrients
              ,Detritus = NULL # no function here
              ,Import   = NULL # from change in MLD and sea ice
              ,Export   = NULL # from change in MLD and sea ice and scavenging
              ) # end pDi 
  ,pSm = list( Name  = "Phyto-small"
                  ,X0     = 20 # possible starting value 
                  ,Attr   = list(Units = "g.m-2" 
                               ,MassToMole = 1/12  # ensure correct units 1 g : 12 mole (if kg then = 1/12000 etc)
                               ,WtC   = NULL     # g.Carbon individual
                               ,WW_C  = NULL     # g.WetWeight (g.carbon)-1
                               ,r_FeC = 0.005   # micromole Fe (millimole C)-1 Hauck
                               ) # end Attributes
               ,Consume = # for phytoplankton, consumption requires an estimate of production (the growth function below simply requires a conversion to carbon)
                 list(params = list(food   = c("nFeM") # names of taxa being consumed - from names(a)
                                    ,MuMax = 0.66 # from Jeffery
                                    ,alpha = 0.088    # (W m-2 d)-1 Hauck - note alpha*q{Chl:N}*q{N:C}
                                    ,k    =  c(0.02)    # millimole Fe m-3 Hauck  - vector of half saturation constants for each food
                                    ) # end list
                    ,fn     = "fnPhConsume")
               ,Produce  = NULL # no production of nutrients
               ,Detritus = NULL # no function here
               ,Import   = NULL # from change in MLD and sea ice
               ,Export   = NULL # from change in MLD and sea ice and scavenging
  ) # end pSm
       ) # end a

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
load(paste0(fPathIn,"Fe_means.Rdata")) # df Iron concentration (micromol.m-3) rows(MEASO area) cols(n, name, value)
load(paste0(fPathIn,"Si_means.Rdata")) # df Silicate concentration (mmol.m-3) rows(MEASO area) cols(n, name, value)
load(paste0(fPathIn,"Lat_means.Rdata")) # df Latitude mean for MEASO areas rows(MEASO area) cols(n, name, value)
load(paste0(fPathIn,"Insol_means.Rdata"))  # matrix Insolation (W.m-2.d-1) rows(months) cols(MEASO area)
load(paste0(fPathIn,"MLD_means.Rdata"))  # matrix Mixed Layer Depth (m) rows(months) cols(MEASO area)
load(paste0(fPathIn,"CICE_means.Rdata"))  # matrix Sea ice percent cover (%) rows(months) cols(MEASO area) 
load(paste0(fPathIn,"Temp_MLD_means.Rdata"))  # matrix Temperature mixed layer (oC) rows(months) cols(MEASO area) 
load(paste0(fPathIn,"Temp_Deep_means.Rdata"))  # matrix Temperature MLD to 1000m (oC) rows(months) cols(MEASO area) 


# 2.1 Managing time ####
DaysInMonth<-c(31,28,31,30,31,30,31,31,30,31,30,31)
DayMidMonth<-c(15,14,15,15,15,15,15,15,15,15,15,15)

