# SOwebFns.r

# by Andrew Constable
# December 2023
# last update 20240206


# Primary Production ####

#     Evans-Parslow routine modified from Melbourne-Thomas et al 2015 ####
#          Angle of incident radiation is dependent on time of day and day of the year (tilt of the Earth) and influences the 
#          incident radiation reaching depths.  Corrections for this are included in the description on page 66 in Eq 5.
 
EvansParslow <- function( # JMT Equation 5 Average Growth Rate given latitude, day length and mixed layer depth
                          # note Jmax(T) from Jeffery et al
   Ph     # character name for phytoplankton group in list
  ,a      # biological parameters
  ,cE     # environmental constants
  ,tE     # environmental time series
  ,PAR    # while PAR is present as a list if cE, only one is parsed.  The one parsed is either with or without sea ice present.
  ,Jmax   # list of timestep vectors of Jmax for different phytoplankton groups  
                 ) { # start function
#  using the following data:
#   tE$day    # Julian day of the year
#  ,tE$Insol    # solar radiation at surface (daily forcing variable)
#  ,tE$MLD    # mixed layer depth (daily forcing variable)
#  ,a$Ph[[Ph]]$alpha  # photosynthetic efficiency
#  ,cE$kw     # constant (Table 1) - light attenuation in water   
#  ,Jmax[[Ph]]   # maximum growth rate as 'nutrient limited (light saturated) growth'
#  ,cE$Lat    # latitude
  ##    
  ## Compute daily averaged light-limited growth rate 
  ## analytical integration over layer thickness and day after Evans and
  ## Parslow (1985)
  #fx1 <- 2*pi*(day+192)/365
  fx1 <- 2*pi*tE$day/365
  declin <- 0.006918-  
    (0.399912*cos(fx1)+0.006758*cos(2*fx1)+0.002697*cos(3*fx1)) +
      +(0.070257*sin(fx1)+0.000907*sin(2*fx1)+0.001480*sin(3*fx1))

  ## Compute solar angle at noon (and assume that this is the
  ## equivalent daily averaged incidence angle for direct+diffuse
  ## radiation). cobeta is cos(incidence angle of solar radiation at
  ## noon)
  fx1 <- pi/180*cE$Lat
  cobeta <- pmax(sin(fx1)*sin(declin)+cos(fx1)*cos(declin),0)
  cobeta <- sqrt(1-(1-cobeta^2)/1.33^2)

  ## Length of day
  fx2 <- pmax(pmin(-tan(fx1)*tan(declin),1),-1)
  daylen <- pmax(1.0e-12,acos(fx2)/pi)

  rayb <- pmax(1.0e-15,exp(-tE$MLD/(cobeta*cE$kw)))

  daylen <- daylen/2
  radbio <- pmax(1.0,PAR*tE$Insol)
  vpbio <- Jmax[[Ph]]
  fx1 <- daylen^2*vpbio/(a$Ph[[Ph]]$alpha*radbio)
  fx3 <- fx1
  fx4 <- fx1/rayb
  fu1 <- sqrt(fx3^2+daylen^2)
  fu1 <- fu1-daylen*log((daylen+fu1)/fx3)
  fu2 <- sqrt(fx4^2+daylen^2)
  fu2 <- fu2-daylen*log((daylen+fu2)/fx4)
  J <- -2*vpbio/tE$MLD*cobeta*cE$kw*(fu1-fu2-fx3+fx4)  # note cE$kw is 1/attenuation coefficient

  J
}
