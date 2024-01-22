## Integrate NPZD system by Euler stepping
NPZD.sim <- function(a,x1,u,day,sra,mld,no3e,w,sst,lat,dt) {
  n <- length(u)+1
  X <- matrix(0,4,n)

  Jmax <- a[3]*a[4]^sst
  J <- EvansParslow(day,sra,mld,u,a[1],a[2],Jmax,-lat)

  X[,1] <- x1
  for(k in 2:n) {
    JA <- Jmax[k-1]*X[1,k-1]/(a[5]+X[1,k-1])
    JB <- min(J[k-1],JA)
    X[,k] <- X[,k-1] + dt*F(k-1,X[,k-1],u[k-1],a,mld[k-1],no3e[k-1],w[k-1],JB)
  }

  X
}


