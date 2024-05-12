# CCAMLR Decision Rules - Input data

ProjYears<-list(YearsTotal = 60
                ,Fishery    = c(11:30) # accounting for Year 0 in dataframe
                ,Recovery   = 20
                ,GammaCalcs = 20
) # end ProjYears

# 1.0 Environment ####

Env<-list( name = "Env"
           ,plot = list(col = "#993300")
           ,Period       = list(P1 = list( years = 10
                                           ,doFn = FALSE
                                           ,fn = NULL
                                           ,data = 100)
                                ,P2 = list( years = 10
                                            ,doFn = FALSE
                                            ,fn = NULL
                                            ,data = 100)
                                ,P3 = list( years = 20
                                            ,doFn = TRUE
                                            ,fn = function(yrs,d){t<-c(1:yrs);100-(t-1)*1.0}
                                            ,data = 100)
                                ,P4 = list( years = 60
                                            ,doFn = TRUE
                                            ,fn = function(yrs,d){rep(80,yrs)}
                                            ,data = 100)
           ) # end Period data
) # end environment

# 2.0 Krill ####

Krill<-list( name="Krill"
             ,plot = list(col = "#FF9900")
             ,par = list(K = Env$Period$P1$data
                       , r=3.1     # 3
                       , phi=1 # 0.1
                       , M = NA
                       , useMrate=FALSE
                       , Mmax=NA
                       , Holling = list(p50 = 0.3   # 0.3   prey at which 50% consumption without predator competition, relative to prey population size at initial PBmax
                                        ,q = 2   )    # 3.5   Holling curvature
                       ) # end par
            ,B0asPropKrillMaxProd=0.766 # 0.97
            ,ePB = 3.09 # ecopath Production to Biomass ratio: Hill et al 2021 Euphausids
            ,modelK = list(changingK=TRUE)
) # end Krill

# 3.0 Predators ####
Predators<-list(
  Whale = list( name="Whale"
            ,plot = list(col = "#0033FF")
            ,QBhat=NA  # calculate from prey maximum production
            ,useHolling = TRUE
            ,Holling = list(p50 = 0.2   # 0.3   prey at which 50% consumption without predator competition, relative to prey population size at initial PBmax
                            ,q = 2       # 3.5   Holling curvature
                            ,c=1         # relative impacts of predators on availability of prey
                            ,g=3         # g>=1 delay in decline of effect of predator competition 
                            ,doPC=FALSE) # adjust Holling for predator competition
            ,BrelP = 0.0318  # biomass relative to krill - McCormack 
            ,PB    = 0.03 # baleen whales from Hill et al 2021
            ,QB    = 5.31 # baleen whales from Hill et al 2021
            ,RB    = NA # baleen whales from Hill et al 2021
            ,A     = 0.84   # McCormack baleen whales from Hill et al 2021
             ) # end whale 

 ,Penguin = list(name="Penguin"
              ,plot = list(col = "#336633")
              ,QBhat=NA  # calculate from prey maximum production
             ,useHolling = TRUE
             ,Holling = list(p50 = 0.3   # 0.3   prey at which 50% consumption without predator competition, relative to prey population size at initial PBmax
                            ,q = 2       # 3.5   Holling curvature
                            ,c=1         # relative impacts of predators on availability of prey
                            ,g=3         # g>=1 delay in decline of effect of predator competition 
                            ,doPC=FALSE) # adjust Holling for predator competition
             ,BrelP = 0.004  # biomass relative to krill - McCormack 
             ,PB    = 0.19 # from Hill et al 2021
             ,QB    = 30.1 # from Hill et al 2021
             ,RB    = NA # from Hill et al 2021
             ,A     = 0.77   # from Hill et al 2021
             ) # end penguin 
 ) # end predators

# 4.0 Fishery ####

Fishery <- list( name="Fishery"
                 ,plot = list(col = "#CC0033")
                ,K_target_0        = 0.75
                ,K_critical_0      = 0.2
                ,Time_for_recovery = 25
) # end fishery

