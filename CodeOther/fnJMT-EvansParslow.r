## NPZD optimal control
## Allows alpha as a control

# edits/revisions by Andrew Constable (AC) January 2024
#   * changed EvansParlsow to EvansParslow
#   * commenting the inputs to EvansParslow routine


EvansParslow <- function( # Equation 5 Average Growth Rate given latitude, day length and mixed layer depth
                  day    # Julian day of the year
                 ,sra    # solar radiation at surface (daily forcing variable)
                 ,mld    # mixed layer depth (daily forcing variable)
                 ,alpha  # photosynthetic efficiency
                 ,par    # PAR as proportion of incident solar radiation immediately below the surface
                 ,kw     # constant (Table 1) - light attenuation in water   
                 ,Jmax   # maximum growth rate as 'nutrient limited (light saturated) growth'
                 ,lat    # latitude
                 ) { # start function

  ## Angle of incident radiation is dependent on time of day and day of the year (tilt of the Earth) and influences the 
  ## incident radiation reaching depths.  Corrections for this are included in the description on page 66 in Eq 5.
  
  ##    
  ## Compute daily averaged light-limited growth rate 
  ## analytical integration over layer thickness and day after Evans and
  ## Parslow (1985)
  #fx1 <- 2*pi*(day+192)/365
  fx1 <- 2*pi*day/365
  declin <- 0.006918-  
    (0.399912*cos(fx1)+0.006758*cos(2*fx1)+0.002697*cos(3*fx1)) +
      +(0.070257*sin(fx1)+0.000907*sin(2*fx1)+0.001480*sin(3*fx1))

  ## Compute solar angle at noon (and assume that this is the
  ## equivalent daily averaged incidence angle for direct+diffuse
  ## radiation). cobeta is cos(incidence angle of solar radiation at
  ## noon)
  fx1 <- pi/180*lat
  cobeta <- pmax(sin(fx1)*sin(declin)+cos(fx1)*cos(declin),0)
  cobeta <- sqrt(1-(1-cobeta^2)/1.33^2)

  ## Length of day
  fx2 <- pmax(pmin(-tan(fx1)*tan(declin),1),-1)
  daylen <- pmax(1.0e-12,acos(fx2)/pi)

  rayb <- pmax(1.0e-15,exp(-mld/(cobeta*kw)))

  daylen <- daylen/2
  radbio <- pmax(1.0,par*sra)
  vpbio <- Jmax
  fx1 <- daylen^2*vpbio/(alpha*radbio)
  fx3 <- fx1
  fx4 <- fx1/rayb
  fu1 <- sqrt(fx3^2+daylen^2)
  fu1 <- fu1-daylen*log((daylen+fu1)/fx3)
  fu2 <- sqrt(fx4^2+daylen^2)
  fu2 <- fu2-daylen*log((daylen+fu2)/fx4)
  J <- -2*vpbio/mld*cobeta*kw*(fu1-fu2-fx3+fx4)

  J
}


