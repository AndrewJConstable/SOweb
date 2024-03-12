# SOwebData.r
# Input data and parameters

a<-list( # list of taxa with their parameters and functions as needed (first letter of name is group, next 2 letters are taxon, 4th & more letter relates to subpool)

# Nutrients ####         
  nFeM =  list(Name    = "Nut-Iron-ML"
              ,X0       = 0     #  starting value
              ,Attr    = list( Units = "umole.m-2" # micromole
                               ,Layer = "M"
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
               ,X0       = 0 # possible starting value
               ,Attr    = list( Units = "umole.m-2" # micromole
                                ,Layer = "D"
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
  ,nFeSI =  list(Name     = "Nut-Iron-SeaIce"
                 ,X0      = 0 # possible starting value
                 ,Attr    = list(Units="umol.m-2"
                                 ,Layer = "SI"
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
               ,X0      = 0 # possible starting value
               ,Attr    = list( Units = "mmole.m-2" # millimole
                                ,Layer = "M"
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
               ,X0       = 0 # possible starting value
               ,Attr    = list( Units = "mmole.m-2" # millimole
                                ,Layer = "D"
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
               ,X0       = 0 # possible starting value
               ,Attr    = list( Units = "mmole.m-2" # millimole
                                ,Layer = "SI"
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

# Detrital carbon ####
 ,dCaSI =  list(Name    = "Detritus-Carbon-Sea ice"
                ,X0       = 0 # possible starting value
                ,Attr    = list( Units = "g.m-2" # grams : 1.2E4 mmole.m-2
                                 ,Layer = "SI"
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
               ,X0       = 0 # possible starting value
               ,Attr    = list( Units = "g.m-2" # grams : 1.2E4 mmole.m-2
                                ,Layer = "M"
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
              ,X0       = 0 # possible starting value
              ,Attr    = list( Units = "g.m-2" # grams : 1.2E4 mmole.m-2
                               ,Layer = "D"
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
              ,X0       = 0 # possible starting value
              ,Attr    = list( Units = "g.m-2" # grams : 1.2E4 mmole.m-2
                               ,Layer = "S"
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


# Detrital Iron ####
# - note that iron from detritus emanating from plants and animals is produced from detrital carbon consumption.

,dFeSI =  list(Name    = "Detritus-Iron-Sea ice"
                ,X0       = 0 # possible starting value
                ,Attr    = list( Units = "umol.m-2" # micromole
                                 ,Layer = "SI"
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
               ,X0       = 0 # possible starting value
               ,Attr    = list( Units = "umol.m-2" # micromole
                                ,Layer = "M"
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
              ,X0       = 0 # possible starting value
              ,Attr    = list( Units = "umol.m-2" # micromole
                               ,Layer = "D"
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
              ,X0       = 0 # possible starting value
              ,Attr    = list( Units = "umol.m-2" # micromole
                               ,Layer = "S"
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

# Detrital Silica ####
# - note that Si from detritus emanating from plants and animals is produced from detrital carbon consumption.

,dSiSI =  list(Name    = "Detritus-Silica-Sea ice"
               ,X0       = 0 # possible starting value
               ,Attr    = list( Units = "mmol.m-2" # millimole
                                ,Layer = "SI"
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
              ,X0       = 0 # possible starting value
              ,Attr    = list( Units = "mmol.m-2" # millimole
                               ,Layer = "M"
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
              ,X0       = 0 # possible starting value
              ,Attr    = list( Units = "mmol.m-2" # millimole
                               ,Layer = "D"
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
              ,X0       = 0 # possible starting value
              ,Attr    = list( Units = "mmol.m-2" # millimole
                               ,Layer = "S"
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


# Phytoplankton ####

  ,pDi = list( Name    = "Phyto-diatom" 
              ,X0       = 100 # possible starting value
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
                  ,X0     = 10 # possible starting value 
                  ,Attr   = list(Units = "g.m-2" 
                               ,WtC   = NULL     # g.Carbon individual
                               ,WW_C  = NULL     # g.WetWeight (g.carbon)-1
                               ,r_FeC = 0.005   # micromole Fe (millimole C)-1 Hauck
                               ,PBratio = 3.0
                               ) # end Attributes
               ,Consume = # for phytoplankton, consumption requires an estimate of production (the growth function below simply requires a conversion to carbon)
                 list(params = list(pool   = c("nFeM") # names of taxa being consumed - from names(a)
                                    ,MuMax = 0.63 # from Jeffery
                                    ,alpha = 0.088    # (W m-2 d)-1 Hauck - note alpha*q{Chl:N}*q{N:C}
                                    ,k    =  c(0.02)    # umole Fe m-3 Hauck  - vector of half saturation constants for each pool
                                    ) # end list
                    ,fn     = "fnPhConsume")
               ,Produce  =  # convert the amount consumed into carbon mass - based on iron uptake
                              list(params = NULL # end list
                                  ,fn     = "fnPhProduceFe")
  ) # end pSm

# Zooplankton ####
,zMi = list( Name  = "Zooplankton-micro"
             ,X0     = 1 # possible starting value 
             ,Attr   = list(Units = "g.m-2" 
                            ,WtC   = NULL     # g.Carbon individual
                            ,WW_C  = NULL     # g.WetWeight (g.carbon)-1
                            ,r_FeC = 0.005   # micromole Fe (millimole C)-1 Hauck
                            ,PBratio = 2.0
                            ,Foraging = list(params = c(0,-300) # Depth range - min, max for determining overlap with consumers
                                             ,fn    =  "fnH_forageUniform")
             ) # end Attributes
             ,Consume = # 
               list(params = list(pool   = c("pSm") # names of taxa being consumed - from names(a)
                                  ,foodSel = {res<-c(1); names(res)<-c("pSm");res} # selectivity of each resource type if encountered
                                  ,actions = list(  pSm  = list(params = list( I    = NULL   # maximum ingestion rate 
                                                                              ,s = NULL)  # selectivity - probability of consuming resource given encounter
                                                                 ,fn     = "fnH_consumeHolling2_KY")
                                  )  # end actions
               ) # end params list
               ,fn     = "fnD_Consume")
             ,Produce  =  # transfer consumption (for detritus need to use mole to mole and carbon to mole functions)
               list(params = list(pool   = c("dSiM","dSiS","pDi") # names of taxa being consumed - from names(a)
                                  ,actions = list(pSm = list(params = list( A    = 1 # assimilation rate: (1-A) goes to detritus 
                                                                           ,conv = 1 # conversion efficiency with remainder not going to detritus
                                                                           ) # end params
                                                              ,fn     = "fnH_produceB")
                                  ) # end actions list
               ) # end params list
               ,fn     = "fnD_Produce") # end produce
) # end zMi

# Loss pools ####
,lRes = list(Name    = "Respiration"
              ,X0       = 0 # possible starting value
              ,Attr    = list( Units = "g.m-2"
              ) # end attributes 
              ,Consume = # 
                list(params = list(pool   = c("pDi","pSm","zMi") # names of taxa being consumed - from names(a)
                                   ,actions = list( pDi  = list(params = list(respCoeff = NULL )  
                                                                ,fn     = "fnGeneralRespire")
                                                    ,pSm = list(params = list(respCoeff = NULL )  
                                                                ,fn     = "fnGeneralRespire")
                                                    ,zMi = list(params = list(respCoeff = NULL )  
                                                                ,fn     = "fnGeneralRespire")
                                   ) # end actions
                ) # end params list
                ,fn     = "fnL_Consume")
) # end lRes  
,lMort = list(Name    = "Mortality not included in production"
              ,X0       = 0 # possible starting value
              ,Attr    = list( Units = "d-1"
                             ) # end attributes 
              ,Consume = # 
                list(params = list(pool   = c("pDi","pSm") # names of taxa being consumed - from names(a)
                                   ,actions = list( pDi  = list(params = list(Rmort = 0.99/200 )  # d-1 
                                                                 ,fn     = "fnGeneralMortality")
                                                    ,pSm = list(params = list(Rmort = 0.99/50 )  # d-1 
                                                                ,fn     = "fnGeneralMortality")
                                   ) # end actions
                ) # end params list
                ,fn     = "fnL_Consume")
             ) # end lMort  

# end A list ####
) # end a

