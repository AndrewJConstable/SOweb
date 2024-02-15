# SOwebData_cE.r
# Input data and parameters
#    1.3. Constants for the environment

cE<-list( par = list(w = 0.43, si = 0.1) # proportion of insolation available as PAR without (w) and with sea ice present (si) sea ice has attenuation of 1.5 m-1
         ,kw  = 1/0.04                     # light attenuation (note that this is converted to 1/0.04 for computational efficiency)
         ,DepthMax = 2000                # maximum depth of deep stratum
         # carbon concentrations are modelled as millimole following Hauck et al 2013.  Thus the conversion relates gC:mmoleC
         # where 12 mole Carbon = 1 gram Carbon
         ,C_MoleToMass = 1/1.2E4            #  converting millimole to grams
         ,C_MassToMole = 1.2E4              # converting gram to millimole
         ,Fe_scavenge = list(Rate = 0.0156*0.01    # ((mmol C m-3)-1 d-1) (Hauck et al 2013) x proportion of DFe which is free Fe (Smith et al 2022)
                            ,FreePropDFe = 0.01    # (Smith et al 2022)
                                      ) # end Fe_scavenge
         ) # end eC

