# SOwebData.r
# Input data and parameters
#    1.1. Minimisation inputs ####
# time0 relates to 1st element of a temporal vector.  The matrix of state variables 
nYrs<-1
nStepsPerYr<-365 
Date1<-as.Date("2015-02-14")  # start date chosen to be at minimum of sea ice to enable incorporation of phytoplankton in sea ice etc.


#    1.2. Input Parameters ####

a<-list( # list of taxa with their parameters and functions as needed (first letter of name is group, next 2 letters are taxon, 4th & more letter relates to subpool)
       
  nFeM =  list(Name    = "Nut-Iron-ML"
              ,X0       = 100 # possible starting value
              ,Attr    = list( Units = "umole.m-2" # micromole
                              ,Which_C_ratio = "r_FeC"
                          )  
              ,Consume = # 
                list(params = list(pool   = c("nFeD","nFeSI","dFeM") # names of taxa being consumed - from names(a)
                                   ,actions = list(nFeD  = list(params = list(do_dMLDgt0 = TRUE)
                                                              ,fn     = "fnN_changeMLD") 
                                                 ,nFeSI = list(params = list(do_dCICEgt0 = FALSE)
                                                              ,fn     = "fnN_changeSI")
                                                 ,dFeM  = list(params = NULL # Hauck 2013 see equations
                                                             ,fn     = "fnN_remineraliseD")
                                                   ) # end actions list
                                  ) # end params list
                      ,fn     = "fnN_Consume")
              ,Produce  =  # transfer consumption
                          list(params = NULL # end list
                              ,fn     = "fnN_Produce")
  ) # end nFeM
 ,nFeD =  list(Name    = "Nut-Iron-Deep"
               ,X0       = 2000 # possible starting value
               ,Attr    = list( Units = "umole.m-2" # micromole
                                ,Which_C_ratio = "r_FeC"
               )  
               ,Consume = # 
                 list(params = list(pool   = c("nFeM","dFeD") # names of taxa being consumed - from names(a)
                                    ,actions = list(nFeM   = list(params = list(do_dMLDgt0 = FALSE)
                                                                 ,fn     = "fnN_changeMLD") 
                                                    ,dFeD  = list(params = NULL # Hauck 2013 see equations
                                                                   ,fn     = "fnN_remineraliseD")
                                    ) # end actions list
                 ) # end params list
                 ,fn     = "fnN_Consume")
               ,Produce  =  # transfer consumption
                 list(params = NULL # end list
                      ,fn     = "fnN_Produce")
  ) # end nFeD
  ,nFeSI =  list(Name    = "Nut-Iron-SeaIce"
                 ,X0       = 100 # possible starting value
                 ,Attr    = list(Units="umol.m-2"
                                 ,Which_C_ratio = "r_FeC"
                 )        
                 ,Consume = # 
                   list(params = list(pool   = c("nFeM") # names of taxa being consumed - from names(a)
                                      ,actions = list(nFeM = list(params = list(do_dCICEgt0 = TRUE)
                                                                    ,fn     = "fnN_changeSI")
                                      ) # end actions list
                   ) # end params list
                   ,fn     = "fnN_Consume")
                 ,Produce  =  # transfer consumption
                   list(params = NULL # end list
                        ,fn     = "fnN_Produce")
  ) # end nFeSI

 
  ,nSiM =  list(Name    = "Nut-Silicate-MLD"
               ,X0       = 4000 # possible starting value
               ,Attr    = list( Units = "mmole.m-2" # millimole
                                ,Which_C_ratio = "r_SiC"
               )  
               ,Consume = # from detrital pool, deep water, sea ice
                 list(params = list(pool   = c("nSiD","nSiSI","dSiM") # names of taxa being consumed - from names(a)
                                  ,actions = list(nSiD  = list(params = list(do_dMLDgt0 = TRUE)
                                                               ,fn     = "fnN_changeMLD") 
                                                  ,nSiSI = list(params = list(do_dCICEgt0 = FALSE)
                                                                ,fn     = "fnN_changeSI")
                                                  ,dSiM  = list(params = NULL # Hauck 2013 C degradation (0.15 d-1) x C reminerisation (0.1 d-1) rates
                                                                ,fn     = "fnN_remineraliseD")
                                  ) # end actions list
                               ) # end params list
                    ,fn     = "fnN_Consume")
               ,Produce  =  # transfer consumption
                 list(params = NULL # end list
                      ,fn     = "fnN_Produce")
  ) # end nSiM
 ,nSiD =  list(Name    = "Nut-Silicate-Deep"
               ,X0       = 80000 # possible starting value
               ,Attr    = list( Units = "mmole.m-2" # millimole
                                ,Which_C_ratio = "r_SiC"
               )  
                 ,Consume = # 
                 list(params = list(pool   = c("nSiM","dSiD") # names of taxa being consumed - from names(a)
                                    ,actions = list(nSiM   = list(params = list(do_dMLDgt0 = FALSE)
                                                                  ,fn     = "fnN_changeMLD") 
                                                    ,dSiD  = list(params = NULL # Hauck 2013 see equations
                                                                  ,fn     = "fnN_remineraliseD")
                                    ) # end actions list
                 ) # end params list
                 ,fn     = "fnN_Consume")
               ,Produce  =  # transfer consumption
                 list(params = NULL # end list
                      ,fn     = "fnN_Produce")
 ) # end nSiD
 ,nSiSI =  list(Name    = "Nut-Silicate-Sea ice"
               ,X0       = 20 # possible starting value
               ,Attr    = list( Units = "mmole.m-2" # millimole
                                ,Which_C_ratio = "r_SiC"
               )  
               ,Consume = # 
                 list(params = list(pool   = c("nSiM") # names of taxa being consumed - from names(a)
                                    ,actions = list(nSiM = list(params = list(do_dCICEgt0 = TRUE)
                                                                ,fn     = "fnN_changeSI")
                                    ) # end actions list
                 ) # end params list
                 ,fn     = "fnN_Consume")
               ,Produce  =  # transfer consumption
                 list(params = NULL # end list
                      ,fn     = "fnN_Produce")
 ) # end nSiSI

 # Detritus ###################################################
 
 # detritus Fe
,dFeSI =  list(Name    = "Detritus-Iron-Sea ice"
                ,X0       = 20 # possible starting value
                ,Attr    = list( Units = "umol.m-2" # micromole
                                 ,Which_C_ratio = "r_FeC"   # for origin of detritus
                               ) # end attributes 
                ,Consume = # 
                  list(params = list(pool   = c("dFeM","pDi","pSm") # names of taxa being consumed - from names(a), repeat name if different functions
                                     ,actions = list(dFeM = list(params = list(do_dCICEgt0 = TRUE)  # umole to umole
                                                                 ,fn     = "fnN_changeSI")
                                                    ,pDi = list(params = list(do_dCICEgt0 = TRUE)   # consume carbon then produce umole
                                                                 ,fn     = "fnN_changeSI")  # end actions list
                                                    ,pSm = list(params = list(do_dCICEgt0 = TRUE)    # consume carbon then produce umole
                                                                 ,fn     = "fnN_changeSI") # end actions list
                                                       ) # end params list
                  ,fn     = "fnD_Consume")
                ,Produce  =  # transfer consumption (for detritus need to use mole to mole and carbon to mole functions)
                  list(params = list(pool   = c("dFeM","pDi","pSm") # names of taxa being consumed - from names(a), repeat name if different functions
                                    ,actions = list(dFeM = list(params = NULL  # umole to umole
                                                                ,fn     = "fnD_produceSameUnits")
                                                   ,pDi = list(params = NULL   # consume carbon then produce umole
                                                                ,fn     = "fnD_produceMoleFromCarbon") ) # end actions list
                                                   ,pSm = list(params = NULL    # consume carbon then produce umole
                                                                ,fn     = "fnD_produceMoleFromCarbon") ) # end actions list
                                  ) # end params list
                       ,fn     = "fnD_Produce"
                  )
 ) # end dFeSI
 
# Phytoplankton ###################################################

  ,pDi = list( Name    = "Phyto-diatom" 
              ,X0       = 30 # possible starting value
              ,Attr   = list( Units = "g.m-2"
                             ,WtC   = NULL     # g.Carbon individual weight
                             ,WW_C  = NULL     # g.WetWeight per carbon (g.carbon)-1
                             ,r_FeC = 0.005   # micromole Fe (millimole C)-1 Hauck
                             ,r_SiC = 0.8     # mole Si (mole C)-1 Hauck
                             ) # end Attributes
              ,Consume = # for phytoplankton, consumption requires an estimate of production (the growth function below simply requires a conversion to carbon)
                              list(params = list(pool   = c("nFeM","nSiM") # names of taxa being consumed - from names(a)
                                                ,MuMax = 1.44 # from Jeffery
                                                ,alpha = 0.16            # (W m-2 d)-1 Hauck - note alpha*q{Chl:N}*q{N:C}
                                                ,k    =  c(0.12,4.0)    # micromole.m-3, millimole.m-3 Hauck  - vector of half saturation constants for each pool
                                                 ) # end list
                                  ,fn     = "fnPhConsume")
              ,Produce  =  # convert the amount consumed into carbon mass - based on iron uptake
                              list(params = NULL # end list
                                  ,fn     = "fnPhProduceFe")
              ) # end pDi 
  ,pSm = list( Name  = "Phyto-small"
                  ,X0     = 20 # possible starting value 
                  ,Attr   = list(Units = "g.m-2" 
                               ,WtC   = NULL     # g.Carbon individual
                               ,WW_C  = NULL     # g.WetWeight (g.carbon)-1
                               ,r_FeC = 0.005   # micromole Fe (millimole C)-1 Hauck
                               ) # end Attributes
               ,Consume = # for phytoplankton, consumption requires an estimate of production (the growth function below simply requires a conversion to carbon)
                 list(params = list(pool   = c("nFeM") # names of taxa being consumed - from names(a)
                                    ,MuMax = 0.66 # from Jeffery
                                    ,alpha = 0.088    # (W m-2 d)-1 Hauck - note alpha*q{Chl:N}*q{N:C}
                                    ,k    =  c(0.02)    # millimole Fe m-3 Hauck  - vector of half saturation constants for each pool
                                    ) # end list
                    ,fn     = "fnPhConsume")
               ,Produce  =  # convert the amount consumed into carbon mass - based on iron uptake
                              list(params = NULL # end list
                                  ,fn     = "fnPhProduceFe")
  ) # end pSm
       ) # end a

#    1.3. Constants for the environment

cE<-list(Lat = eLat
         ,par = list(w = 0.43, si = 0.1) # proportion of insolation available as PAR without (w) and with sea ice present (si) sea ice has attenuation of 1.5 m-1
         ,kw  = 0.04                     # light attenuation (note that this is converted to 1/0.04 for computational efficiency)
         ,Fe  = nFe_0                    # initial iron concentration in deep water
         ,Si  = nSi_0                    # initial silicic acid concentration in deep water
         ,DepthMax = 2000                # maximum depth of deep stratum
         ,C_MoleToMass = 1/12            # Carbon mole to mass in 12 mole to 1 g carbon ((if kg then = 1/12000 etc))
         ,C_MassToMole = 12      
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

