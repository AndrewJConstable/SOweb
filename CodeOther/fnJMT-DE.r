## RHS of NPZD system
F <- function(k # vector element to read
             ,x # X vector of NPZD
             ,u # mortality rate
             ,a # parameters
             ,mld # mixed layer depth for time step
             ,no3e # nitrate concentration for time step
             ,w # change in mixed layer depth for time step
             ,J # phytoplankton growth rate
             ) {

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

