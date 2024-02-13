# SOwebData.r
# Input data and parameters
#    1.1. Minimisation inputs ####
# time0 relates to 1st element of a temporal vector.  The matrix of state variables 
nYrs<-1
nStepsPerYr<-365 
Date1<-as.Date("2015-02-14")  # start date chosen to be at minimum of sea ice to enable incorporation of phytoplankton in sea ice etc.


#    1.2. Input Parameters ####

a<-list( # list of taxa with their parameters and functions as needed (first letter of name is group, next 2 letters are taxon, 4th & more letter relates to subpool)

# Nutrients ####         
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

 # detritus carbon
 ,dCaSI =  list(Name    = "Detritus-Carbon-Sea ice"
                ,X0       = 20 # possible starting value
                ,Attr    = list( Units = "g.m-2" # grams
                ) # end attributes 
                ,Consume = # 
                  list(params = list(pool   = c("dCaM","pDi","pSm") # names of taxa being consumed - from names(a)
                                     ,actions = list(dCaM = list(params = list(do_dCICEgt0 = TRUE)  # umole to umole
                                                                 ,fn     = "fnD_changeSI")
                                                     ,pDi = list(params = list(do_dCICEgt0 = TRUE)   # consume carbon then produce umole
                                                                 ,fn     = "fnD_changeSI")  # end actions list
                                                     ,pSm = list(params = list(do_dCICEgt0 = TRUE)    # consume carbon then produce umole
                                                                 ,fn     = "fnD_changeSI") # end actions list
                                     ) # end actions
                                     ) # end params list
                                     ,fn     = "fnD_Consume")
                       ,Produce  =  # transfer consumption (for detritus need to use mole to mole and carbon to mole functions)
                         list(params = list(pool   = c("dCaM","pDi","pSm") # names of taxa being consumed - from names(a), repeat name if different functions
                                            ,actions = list(dCaM = list(params = NULL  # umole to umole
                                                                        ,fn     = "fnD_produceSameUnits")
                                                            ,pDi = list(params = NULL   # consume carbon then produce umole
                                                                        ,fn     = "fnD_produceSameUnits")
                                                            ,pSm = list(params = NULL    # consume carbon then produce umole
                                                                        ,fn     = "fnD_produceSameUnits")
                                            ) # end actions
                         ) # end params list
                       ,fn     = "fnD_Produce") # end produce
 ) # end dCaSI
 
 ,dCaM =  list(Name    = "Detritus-Carbon- ML"
               ,X0       = 20 # possible starting value
               ,Attr    = list( Units = "g.m-2" # grams
                               ) # end attributes 
               ,Consume = # 
                 list(params = list(pool   = c("dCaSI","dCaD","pDi","pSm") # names of taxa being consumed - from names(a)
                                    ,actions = list( dCaSI = list(params = list(do_dCICEgt0 = FALSE)  # umole to umole
                                                                   ,fn     = "fnD_changeSI")
                                                     ,dCaD  = list(params = list(do_dMLDgt0 = TRUE)  # umole to umole
                                                                   ,fn     = "fnD_DeepToML")
                                                     ,pDi = list(params = list(Ragg=0.015)   # aggregation rate Hauck et al 2013
                                                                   ,fn     = "fnD_consumePhMortality")  # end actions list
                                                     ,pSm = list(params = list(Ragg=0.015)   # aggregation rate Hauck et al 2013
                                                                   ,fn     = "fnD_consumePhMortality") # end actions list
                                    ) # end actions
                 ) # end params list
               ,fn     = "fnD_Consume")
 ,Produce  =  # transfer consumption (for detritus need to use mole to mole and carbon to mole functions)
   list(params = list(pool   = c("dCaSI","dCaD","pDi","pSm") # names of taxa being consumed - from names(a)
                      ,actions = list( dCaSI = list(params = NULL  # umole to umole
                                                    ,fn     = "fnD_produceSameUnits")
                                       ,dCaD  = list(params = NULL  # umole to umole
                                                     ,fn     = "fnD_produceSameUnits")
                                       ,pDi = list(params = NULL   
                                                   ,fn     = "fnD_produceSameUnits")  
                                       ,pSm = list(params = NULL    
                                                   ,fn     = "fnD_produceSameUnits") 
                      ) # end actions
                      ) # end params list
                ,fn     = "fnD_Produce") # end produce
) # end dCaM

,dCaD =  list(Name    = "Detritus-Carbon- Deep"
              ,X0       = 20 # possible starting value
              ,Attr    = list( Units = "g.m-2" # grams
                             ) # end attributes 
              ,Consume = # 
                list(params = list(pool   = c("dCaM","dCaS","pDi","pSm") # names of taxa being consumed - from names(a)
                                   ,actions = list( dCaM  = list(params = list(rSink = 5.0)  # m.d-1
                                                                 ,fn     = "fnD_MLtoDeep")
                                                    ,dCaS  = list(params = list(rSuspend = 0.0 )  # m d-1 
                                                                  ,fn     = "fnD_resuspend")
                                                    ,pDi = list(params = list(do_dMLDgt0 = FALSE)  # umole to umole
                                                                ,fn     = "fnN_changeMLD")
                                                    ,pSm = list(params = list(do_dMLDgt0 = FALSE)  # umole to umole
                                                                ,fn     = "fnN_changeMLD")
                                                    ) # end actions
                                   ) # end params list
                                   ,fn     = "fnD_Consume")
                     ,Produce  =  # transfer consumption (for detritus need to use mole to mole and carbon to mole functions)
                       list(params = list(pool   = c("dCaM","dCaS","pDi","pSm") # names of taxa being consumed - from names(a)
                                          ,actions = list( dCaM = list(params = NULL  # umole to umole
                                                                        ,fn     = "fnD_produceSameUnits")
                                                           ,dCaS  = list(params = NULL  # umole to umole
                                                                         ,fn     = "fnD_produceSameUnits")
                                                           ,pDi = list(params = NULL   
                                                                       ,fn     = "fnD_produceSameUnits")  
                                                           ,pSm = list(params = NULL    
                                                                       ,fn     = "fnD_produceSameUnits") 
                                          ) # end actions
                                          ) # end params list
                            ,fn     = "fnD_Produce") # end produce
                       ) # end dCaD
                     
,dCaS =  list(Name    = "Detritus-Carbon- Sediment"
              ,X0       = 20 # possible starting value
              ,Attr    = list( Units = "g.m-2" # grams
              ) # end attributes 
              ,Consume = # 
                list(params = list(pool   = c("dCaD") # names of taxa being consumed - from names(a), repeat name if different functions
                                   ,actions = list( dCaD  = list(params = list(rSink = 5.0 )  # m d-1 
                                                                  ,fn     = "fnD_sinkFromDeep")
                                                  ) # end actions
                                 ) # end params list
                ,fn     = "fnD_Consume")            
              ,Produce  =  # transfer consumption (for detritus need to use mole to mole and carbon to mole functions)
                list(params = list(pool   = c("dCaD") # names of taxa being consumed - from names(a), repeat name if different functions
                                   ,actions = list( dCaD  = list(params = NULL  # umole to umole
                                                                  ,fn     = "fnD_produceSameUnits")
                                                 ) # end actions
                               ) # end params list
                ,fn     = "fnD_Produce") # end produce
) # end dCaD


# Iron - note that iron from detritus emanating from plants and animals is produced from detrital carbon consumption.


  
 # detritus Fe
,dFeSI =  list(Name    = "Detritus-Iron-Sea ice"
                ,X0       = 20 # possible starting value
                ,Attr    = list( Units = "umol.m-2" # micromole
                                 ,Which_C_ratio = "r_FeC"   # for origin of detritus
                               ) # end attributes 
                ,Consume = #  - note consume only includes direct transfer from pools of same type.  otherwise relying on carbon detrital pools for produciton
                  list(params = list(pool   = c("dFeM") # names of taxa being consumed - from names(a)
                                     ,actions = list(dFeM = list(params = list(do_dCICEgt0 = TRUE)  # umole to umole
                                                                 ,fn     = "fnD_changeSI")
                                                       ) # end actions
                                     ) # end params list
                  ,fn     = "fnD_Consume") # end consume
                ,Produce  =  # transfer consumption (for detritus need to use mole to mole and carbon to mole functions)
                  list(params = list(pool   = c("dFeM","pDi","pSm") # names of taxa being consumed - from names(a)
                                    ,actions = list(dFeM = list(params = NULL  # umole to umole
                                                                ,fn     = "fnD_produceSameUnits")
                                                   ,pDi = list(params = list(Cpool = "dCaSI")   # consume carbon then produce umole
                                                                ,fn     = "fnD_produceMoleFromCarbonPool") # end actions list
                                                   ,pSm = list(params = list(Cpool = "dCaSI")    # consume carbon then produce umole
                                                                ,fn     = "fnD_produceMoleFromCarbonPool") 
                                                 ) # end actions list
                                  ) # end params list
                       ,fn     = "fnD_Produce") # end produce
 ) # end dFeSI

,dFeM =  list(Name    = "Detritus-Iron-Mixed Layer"
               ,X0       = 20 # possible starting value
               ,Attr    = list( Units = "umol.m-2" # micromole
                                ,Which_C_ratio = "r_FeC"   # for origin of detritus
               ) # end attributes 
               ,Consume = # 
                 list(params = list(pool   = c("nFeM","dFeSI","dFeD") # names of taxa being consumed - from names(a)
                                    ,actions = list( nFeM  = list(params = list(CarbonScavengingPool = "dCaM")  # scavenging umole to umole
                                                                  ,fn     = "fnD_consumeFeScavenge")
                                                    ,dFeSI = list(params = list(do_dCICEgt0 = FALSE)  # umole to umole
                                                                ,fn     = "fnD_changeSI")
                                                    ,dFeD  = list(params = list(do_dMLDgt0 = TRUE)  # umole to umole
                                                                  ,fn     = "fnD_DeepToML")
                                                ) # end actions
                                    ) # end params list
                                    ,fn     = "fnD_Consume")
                      ,Produce  =  # transfer consumption (for detritus need to use mole to mole and carbon to mole functions)
                        list(params = list(pool   = c("nFeM","dFeSI","dFeD","pDi","pSm") # names of taxa being consumed - from names(a)
                                           ,actions = list(nFeM = list(params = NULL  # umole to umole
                                                                       ,fn     = "fnD_produceSameUnits")
                                                           ,dFeSI = list(params = NULL  # umole to umole
                                                                        ,fn     = "fnD_produceSameUnits")
                                                           ,dFeD = list(params = NULL  # umole to umole
                                                                        ,fn     = "fnD_produceSameUnits")
                                                           ,pDi = list(params = list(Cpool = "dCaM")   # consume carbon then produce umole
                                                                       ,fn     = "fnD_produceMoleFromCarbonPool") # end actions list
                                                           ,pSm = list(params = list(Cpool = "dCaM")    # consume carbon then produce umole
                                                                       ,fn     = "fnD_produceMoleFromCarbonPool") 
                                                           ) # end actions list
                                           ) # end params list
                      ,fn     = "fnD_Produce") # end produce
) # end dFeM

,dFeD =  list(Name    = "Detritus-Iron-Deep"
              ,X0       = 20 # possible starting value
              ,Attr    = list( Units = "umol.m-2" # micromole
                               ,Which_C_ratio = "r_FeC"   # for origin of detritus
              ) # end attributes 
              ,Consume = # 
                list(params = list(pool   = c("nFeD","dFeM","dFeS") # names of taxa being consumed - from names(a)
                                   ,actions = list( nFeD  = list(params = list(CarbonScavengingPool = "dCaD")  # scavenging umole to umole
                                                                 ,fn     = "fnD_consumeFeScavenge")
                                                    ,dFeM  = list(params = list(rSink = 5.0)  # m.d-1
                                                                  ,fn     = "fnD_MLtoDeep")
                                                    ,dFeS  = list(params = list(rSuspend = 0.0 )  # m d-1 
                                                                  ,fn     = "fnD_resuspend")
                                   ) # end actions
                ) # end params list
                ,fn     = "fnD_Consume")
              ,Produce  =  # transfer consumption (for detritus need to use mole to mole and carbon to mole functions)
                list(params = list(pool   = c("nFeD","dFeM","dFeS","pDi","pSm") # names of taxa being consumed - from names(a)
                                   ,actions = list(nFeD = list(params = NULL  # umole to umole
                                                               ,fn     = "fnD_produceSameUnits")
                                                   ,dFeM = list(params = NULL  # umole to umole
                                                                ,fn     = "fnD_produceSameUnits")
                                                   ,dFeS = list(params = NULL  # umole to umole
                                                                ,fn     = "fnD_produceSameUnits")
                                                   ,pDi = list(params = list(Cpool = "dCaD")   # consume carbon then produce umole
                                                               ,fn     = "fnD_produceMoleFromCarbonPool") # end actions list
                                                   ,pSm = list(params = list(Cpool = "dCaD")    # consume carbon then produce umole
                                                               ,fn     = "fnD_produceMoleFromCarbonPool") 
                                                  ) # end actions list
                                  ) # end params list
                ,fn     = "fnD_Produce") # end produce
) # end dFeD

,dFeS =  list(Name    = "Detritus-Iron-Deep"
              ,X0       = 20 # possible starting value
              ,Attr    = list( Units = "umol.m-2" # micromole
                               ,Which_C_ratio = "r_FeC"   # for origin of detritus
              ) # end attributes 
              ,Consume = # 
                list(params = list(pool   = c("dFeD") # names of taxa being consumed - from names(a), repeat name if different functions
                                   ,actions = list( dFeD  = list(params = list(rSink = 5.0 )  # m d-1 
                                                                 ,fn     = "fnD_sinkFromDeep")
                                   ) # end actions
                ) # end params list
                ,fn     = "fnD_Consume")
              ,Produce  =  # transfer consumption (for detritus need to use mole to mole and carbon to mole functions)
                list(params = list(pool   = c("dFeD") # names of taxa being consumed - from names(a), repeat name if different functions
                                   ,actions = list(dFeD = list(params = NULL  # umole to umole
                                                               ,fn     = "fnD_produceSameUnits")
                                   ) # end actions list
                ) # end params list
                ,fn     = "fnD_Produce") # end produce
) # end dFeS

# Silica ##########################################################
# - note that Si from detritus emanating from plants and animals is produced from detrital carbon consumption.

,dSiSI =  list(Name    = "Detritus-Silica-Sea ice"
               ,X0       = 20 # possible starting value
               ,Attr    = list( Units = "umol.m-2" # micromole
                                ,Which_C_ratio = "r_SiC"   # for origin of detritus
               ) # end attributes 
               ,Consume = #  - note consume only includes direct transfer from pools of same type.  otherwise relying on carbon detrital pools for produciton
                 list(params = list(pool   = c("dSiM") # names of taxa being consumed - from names(a)
                                    ,actions = list(dSiM = list(params = list(do_dCICEgt0 = TRUE)  # umole to umole
                                                                ,fn     = "fnD_changeSI")
                                    ) # end actions
                 ) # end params list
                 ,fn     = "fnD_Consume") # end consume
               ,Produce  =  # transfer consumption (for detritus need to use mole to mole and carbon to mole functions)
                 list(params = list(pool   = c("dSiM","pDi") # names of taxa being consumed - from names(a)
                                    ,actions = list(dSiM = list(params = NULL  # umole to umole
                                                                ,fn     = "fnD_produceSameUnits")
                                                    ,pDi = list(params = list(Cpool = "dCaSI")   # consume carbon then produce umole
                                                                ,fn     = "fnD_produceMoleFromCarbonPool") # end actions list
                                    ) # end actions list
                 ) # end params list
                 ,fn     = "fnD_Produce") # end produce
) # end dSiSI

,dSiM =  list(Name    = "Detritus-Silica-Mixed Layer"
              ,X0       = 20 # possible starting value
              ,Attr    = list( Units = "umol.m-2" # micromole
                               ,Which_C_ratio = "r_SiC"   # for origin of detritus
              ) # end attributes 
              ,Consume = # 
                list(params = list(pool   = c("dSiSI","dSiD") # names of taxa being consumed - from names(a)
                                   ,actions = list( dSiSI = list(params = list(do_dCICEgt0 = FALSE)  # umole to umole
                                                                  ,fn     = "fnD_changeSI")
                                                    ,dSiD  = list(params = NULL  # umole to umole
                                                                  ,fn     = "fnD_DeepToML")
                                   ) # end actions
                ) # end params list
                ,fn     = "fnD_Consume")
              ,Produce  =  # transfer consumption (for detritus need to use mole to mole and carbon to mole functions)
                list(params = list(pool   = c("dSiSI","dSiD","pDi") # names of taxa being consumed - from names(a)
                                   ,actions = list( dSiSI = list(params = NULL  # umole to umole
                                                                 ,fn     = "fnD_produceSameUnits")
                                                   ,dSiD = list(params = NULL  # umole to umole
                                                                ,fn     = "fnD_produceSameUnits")
                                                   ,pDi = list(params = list(Cpool = "dCaM")   # consume carbon then produce umole
                                                               ,fn     = "fnD_produceMoleFromCarbonPool") # end actions list
                                   ) # end actions list
                ) # end params list
                ,fn     = "fnD_Produce") # end produce
) # end dSiM

,dSiD =  list(Name    = "Detritus-Silica-Deep"
              ,X0       = 20 # possible starting value
              ,Attr    = list( Units = "umol.m-2" # micromole
                               ,Which_C_ratio = "r_SiC"   # for origin of detritus
              ) # end attributes 
              ,Consume = # 
                list(params = list(pool   = c("dSiM","dSiS") # names of taxa being consumed - from names(a)
                                   ,actions = list(  dSiM  = list(params = list(rSink = 5.0)  # m.d-1
                                                                  ,fn     = "fnD_MLtoDeep")
                                                    ,dSiS  = list(params = list(rSuspend = 0.0 )  # m d-1 
                                                                  ,fn     = "fnD_resuspend")
                                   ) # end actions
                ) # end params list
                ,fn     = "fnD_Consume")
              ,Produce  =  # transfer consumption (for detritus need to use mole to mole and carbon to mole functions)
                list(params = list(pool   = c("dSiM","dSiS","pDi") # names of taxa being consumed - from names(a)
                                   ,actions = list(dSiM = list(params = NULL  # umole to umole
                                                                ,fn     = "fnD_produceSameUnits")
                                                   ,dSiS = list(params = NULL  # umole to umole
                                                                ,fn     = "fnD_produceSameUnits")
                                                   ,pDi = list(params = list(Cpool = "dCaD")   # consume carbon then produce umole
                                                               ,fn     = "fnD_produceMoleFromCarbonPool") # end actions list
                                   ) # end actions list
                ) # end params list
                ,fn     = "fnD_Produce") # end produce
) # end dSiD

,dSiS =  list(Name    = "Detritus-Silica-Deep"
              ,X0       = 20 # possible starting value
              ,Attr    = list( Units = "umol.m-2" # micromole
                               ,Which_C_ratio = "r_SiC"   # for origin of detritus
              ) # end attributes 
              ,Consume = # 
                list(params = list(pool   = c("dSiD") # names of taxa being consumed - from names(a)
                                   ,actions = list( dSiD  = list(params = list(rSink = 5.0 )  # m d-1 
                                                                  ,fn     = "fnD_sinkFromDeep")
                                   ) # end actions
                ) # end params list
                ,fn     = "fnD_Consume")
              ,Produce  =  # transfer consumption (for detritus need to use mole to mole and carbon to mole functions)
                list(params = list(pool   = c("dSiD") # names of taxa being consumed - from names(a)
                                   ,actions = list(dSiD = list(params = NULL  # umole to umole
                                                               ,fn     = "fnD_produceSameUnits")
                                   ) # end actions list
                ) # end params list
                ,fn     = "fnD_Produce") # end produce
) # end dSiS


# Phytoplankton ###################################################

  ,pDi = list( Name    = "Phyto-diatom" 
              ,X0       = 30 # possible starting value
              ,Attr   = list( Units = "g.m-2"
                             ,WtC   = NULL     # g.Carbon individual weight
                             ,WW_C  = NULL     # g.WetWeight per carbon (g.carbon)-1
                             ,r_FeC = 0.005   # micromole Fe (millimole C)-1 Hauck
                             ,r_SiC = 0.8     # mole Si (mole C)-1 Hauck
                             ,PBratio = 3.0
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
                  ,X0     = 300 # possible starting value 
                  ,Attr   = list(Units = "g.m-2" 
                               ,WtC   = NULL     # g.Carbon individual
                               ,WW_C  = NULL     # g.WetWeight (g.carbon)-1
                               ,r_FeC = 0.005   # micromole Fe (millimole C)-1 Hauck
                               ,PBratio = 3.0
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
,zMort = list(Name    = "Mortality not included in production"
              ,X0       = 0 # possible starting value
              ,Attr    = list( Units = "m-2"
                             ) # end attributes 
              ,Consume = # 
                list(params = list(pool   = c("pDi","pSm") # names of taxa being consumed - from names(a)
                                   ,actions = list( pDi  = list(params = list(Rmort = 0.99/365 )  # d-1 
                                                                 ,fn     = "fnGeneralMortality")
                                                    ,pSm = list(params = list(Rmort = 0.99/365 )  # d-1 
                                                                ,fn     = "fnGeneralMortality")
                                   ) # end actions
                ) # end params list
                ,fn     = "fnD_Consume")
             ) # end zMort  
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
         ,Fe_scavenge = list(Rate = 0.0156*0.01    # ((mmol C m-3)-1 d-1) (Hauck et al 2013) x proportion of DFe which is free Fe (Smith et al 2022)
                            ,FreePropDFe = 0.01    # (Smith et al 2022)
                                      ) # end Fe_scavenge
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

