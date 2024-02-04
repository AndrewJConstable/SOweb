Jmax<-1
alpha<-0.6
I<-c(0:1000)/100
J1<-alpha*I/((Jmax^2+(alpha*I)^2)^0.5)
plot(I,J1,type="l")
J2<-1-exp(-alpha*I/Jmax)
lines(I,J2,col="red")

T<--2+c(0,1000)/50
K<-273.15
T1<-exp(-4500*(1/(T+K)-1/288.15))
T2<-0.8*exp(-4000*(1/(T+K)-1/293.15))
T3<-0.6*1.066^T
T4<-0.63*exp(0.06*T)
T5<-1.44*exp(0.06*T)
plot(T,T1,type="l",ylim=c(0,3))
lines(T,T2,col="green")
lines(T,T3,col="red")
lines(T,T4,col="blue")
lines(T,T5,col="orange")

        
H<-c(0,100)/50
H1<-exp(-1.5*H)
H2<-exp(-4*H)
plot(H,H1,type="l")
lines(H,H2,col="green")

dirOut<-"./tif/"
fOut<-"Climatology_Si.tiff"
t1<-rast(paste0(dirOut,fOut))
t1<-crop(t1,ext(-180,180,-90,-60))
t2<-project(t1,"EPSG:3976")

tmp3<-project(tmp2,"EPSG:4326")


