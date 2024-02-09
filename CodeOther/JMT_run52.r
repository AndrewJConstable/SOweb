## This code runs the optimal control algorithm
## to reproduce results shown in Figs. 5 & 6 in
## the main text (to reproduce results in Fig. 6) 
## set latitude <- 58 (line 10)
## 
## Requires Rcgmin package

source("NPZD.r")
library(Rcgmin)
latitude <- 52

## Vector of parameters
a <-
  c(
    ## Phytoplankton coefficients
    par   = 0.43,        #Phyto: Photosynthetically active radiation
    kw    = 1/0.04,      #Phyto: Light attenuation due to water
    a     = 0.6,         #Phyto: Max. growth rate parameter
    b     = 1.066,       #Phyto: Max. growth rate parameter
    K1    = 0.5,         #Phyto: Half-sat. const. for N uptake
    muep  = 0.00,		 #Phyto: Linear mortality
    muepq = 0.03,        #Phyto: Quad. mortality

    ## Zooplankton coefficients
    gam1   = 0.75,       #Zoo: Assim. efficiency
    g      = 1.00,       #Zoo: Max. graz. rate
    eps    = 3.00,       #Zoo: Prey capture rate
    muez   = 0.02,       #Zoo: Quad. mortality
    gam2   = 0.01,       #Zoo: Excretion

    ## Detrital coefficients
    mued   = 0.05,       #Det: Remineralization rate
    w_sink = 5           #Det: Sinking velocity
)

## Read data
d <- read.csv(paste("forcing",latitude,".csv",sep=""),header=T)

################################
## Environmental forcings
Jmax <- a[3]*a[4]^d$SST
J <- EvansParlsow(d$Day,d$SRA,d$MLD,0.025,a[1],a[2],Jmax,-latitude)
opar <- par(mfrow=c(2,2))
plot(MLD~Day,type="l",data=d)
plot(SST~Day,type="l",data=d)
plot(SRA~Day,type="l",data=d)
plot(J~Day,type="l",data=d,ylim=range(c(J,Jmax)))
lines(Jmax~Day,type="l",data=d,col="red")
par(opar)


################################
## Expand to multiple years
nyears <- 3
nstep <- 1
## Repeat periodically
sra <- rep(d$SRA,nyears)
mld <- rep(d$MLD,nyears)
no3e <- rep(d$NO3E,nyears)
sst <- rep(d$SST,nyears)

## Interpolate to sub-day scale
day <- seq(1,nyears*365,by=1/nstep)
sra <- approx(1:length(sra),sra,day)$y
mld <- approx(1:length(mld),mld,day)$y
no3e <- approx(1:length(no3e),no3e,day)$y
sst <- approx(1:length(sst),sst,day)$y
w <- c(nstep*diff(mld),0)


################################
## Observed P
d <- read.csv(paste("observedP",latitude,".csv",sep=""),header=T)

## Expand to multiple years
Pobs <- d$obsP[match(((day-1)%%365+1),d$Day)]


################################
## Guess some starting values
x1 <- c(N=18,P=0.1,Z=0.4,D=0.1)


########################################################################################
## Conjugate gadient minimization

## Starting value for u
fit <- list(par=rep(0.02,length(day)))
pr <- fit$par


################################
## Stage 1
## Increased Z holds P down
smooth <- 1.0E1
cost <- NPZD.costu(a,x1,day,sra,mld,no3e,w,sst,latitude,dt=1/nstep,Pobs,smooth=smooth)
grad <- NPZD.gradu(a,x1,day,sra,mld,no3e,w,sst,latitude,dt=1/nstep,Pobs,smooth=smooth)
fit <- Rcgmin(pr,fn=cost,gr=grad,
              lower=rep(0.00001,length(pr)),
              control=list(trace=1,maxit=100))

pr <- fit$par
Xfit <- t(NPZD.sim(a,x1,pr,day,sra,mld,no3e,w,sst,latitude,dt=1/nstep))
cost(pr)
sum(grad(pr))


################################
## Stage 2

##Update starting values
x1 <- Xfit[2*365+1,]

smooth <- 1.0E1
cost <- NPZD.costu(a,x1,day,sra,mld,no3e,w,sst,latitude,dt=1/nstep,Pobs,smooth=smooth)
grad <- NPZD.gradu(a,x1,day,sra,mld,no3e,w,sst,latitude,dt=1/nstep,Pobs,smooth=smooth)
fit <- Rcgmin(pr,fn=cost,gr=grad,
              lower=rep(0.00001,length(pr)),
              control=list(trace=1,maxit=100))

pr <- fit$par
Xfit <- t(NPZD.sim(a,x1,pr,day,sra,mld,no3e,w,sst,latitude,dt=1/nstep))
cost(pr)
sum(grad(pr))


################################
## Stage 3

##Update starting values
x1 <- Xfit[2*365+1,]

smooth <- 1.0E3
cost <- NPZD.costu(a,x1,day,sra,mld,no3e,w,sst,latitude,dt=1/nstep,Pobs,smooth=smooth)
grad <- NPZD.gradu(a,x1,day,sra,mld,no3e,w,sst,latitude,dt=1/nstep,Pobs,smooth=smooth)
fit <- Rcgmin(pr,fn=cost,gr=grad,
              lower=rep(0.00001,length(pr)),
              control=list(trace=1,maxit=100))

pr <- fit$par
Xfit <- t(NPZD.sim(a,x1,pr,day,sra,mld,no3e,w,sst,latitude,dt=1/nstep))
cost(pr)
sum(grad(pr))


################################
## Stage 4

##Update starting values
x1 <- Xfit[2*365+1,]

## Resmooth
pr <- fit$par
plot(pr,pch=16,cex=0.5)
pr <- fitted(loess(pr~day,span=0.1))
lines(pr,col="red")

smooth <- 1.0E4
cost <- NPZD.costu(a,x1,day,sra,mld,no3e,w,sst,latitude,dt=1/nstep,Pobs,smooth=smooth)
grad <- NPZD.gradu(a,x1,day,sra,mld,no3e,w,sst,latitude,dt=1/nstep,Pobs,smooth=smooth)
fit <- Rcgmin(pr,fn=cost,gr=grad,
              lower=rep(0.00001,length(pr)),
              control=list(trace=1,maxit=100))

pr <- fit$par
Xfit <- t(NPZD.sim(a,x1,pr,day,sra,mld,no3e,w,sst,latitude,dt=1/nstep))
cost(pr)
sum(grad(pr))


################################
## Stage 5

##Update starting values
x1 <- Xfit[2*365+1,]

smooth <- 1.0E6
cost <- NPZD.costu(a,x1,day,sra,mld,no3e,w,sst,latitude,dt=1/nstep,Pobs,smooth=smooth)
grad <- NPZD.gradu(a,x1,day,sra,mld,no3e,w,sst,latitude,dt=1/nstep,Pobs,smooth=smooth)
fit <- Rcgmin(pr,fn=cost,gr=grad,
              lower=rep(0.00001,length(pr)),
              control=list(trace=1,maxit=100))

pr <- fit$par
Xfit <- t(NPZD.sim(a,x1,pr,day,sra,mld,no3e,w,sst,latitude,dt=1/nstep))
cost(pr)
sum(grad(pr))

