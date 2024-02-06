## Integrate NPZD system by Euler stepping
# comments by Andrew Constable

NPZD.sim <- function(
    a       # parameters for functions
   ,x1      # vector of NPZD initial values c(N=18,P=0.1,Z=0.4,D=0.1)
   ,u       # photosynthetic efficiency
   ,day     # Julian day of the year
   ,sra    # solar radiation at surface (daily forcing variable)
   ,mld    # mixed layer depth (daily forcing variable)
   ,no3e   # vector of nitrate concentrations for deep layer
   ,w      # change in mixed layer depth during time step (last element =0)
   ,sst    # vector of SST
   ,lat    # latitude of projection (for incident light angle)
   ,dt     # magnitude of time step (related to each element of vectors)
      ) { # start function
  n <- length(u)+1    # time steps
  X <- matrix(0,4,n)

  Jmax <- a[3]*a[4]^sst # vector of Jmax
  J <- EvansParslow(day,sra,mld,u,a[1],a[2],Jmax,-lat)  # vector of J
  
  X[,1] <- x1
  for(k in 2:n) { # 2:n is based on the matrix of state variables
    JA <- Jmax[k-1]*X[1,k-1]/(a[5]+X[1,k-1])
    JB <- min(J[k-1],JA)
    X[,k] <- X[,k-1] + dt*F(k-1,X[,k-1],u[k-1],a,mld[k-1],no3e[k-1],w[k-1],JB)
  }

  X
}


