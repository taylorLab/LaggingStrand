
###
# Nucleosome sub rate
fr<-c(-125:125)
gidLn<-length(gm$gid)
nucf30o80c30<-nucs$gid[which(nucs$occupancy>80 & nucs$fuzzyness<30 & nucs$tssProximity==0 & (ojrWTrep$forCovDD[nucs$gid]+ojrWTrep$revCovDD[nucs$gid])>30)]


nucSel<-nucf30o80c30
sig<-matrix(c((ojrWTrep$forClobber+ojrWTrep$revClobber)/2,gidLn,1))
nucSelRealSCMTOjr<-signalCounterMatrixTri(nSites=nucSel,focalRange=fr,triTransform=gm$gid,signal=sig,sqDF=gm,mask=gm$exclude)
sig<-matrix(gm$polCall,gidLn,1)
nucSelRealSCMTpolCall<-signalCounterMatrixTri(nSites=nucSel,focalRange=fr,triTransform=gm$gid,signal=sig,sqDF=gm,mask=gm$exclude)
sig<-matrix(gm$polCount,gidLn,1)
nucSelRealSCMTpol<-signalCounterMatrixTri(nSites=nucSel,focalRange=fr,triTransform=gm$gid,signal=sig,sqDF=gm,mask=gm$exclude)


vec2<-c(nucSelRealSCMTpol$NormRate/nucSelRealSCMTpolCall$NormRate)
yMean1<-mean(nucSelRealSCMTOjr$NormRate)
yMean2<-mean(vec2)
yTen<-c(yMean2*.9,yMean2*1.1)
yMarg=.15
yBounds1<-c(yMean1*(1-yMarg),yMean1*(1+yMarg))
yBounds2<-c(yMean2*(1-yMarg),yMean2*(1+yMarg))


plot(nucSelRealSCMTOjr$focal,vec2,ylim=yBounds2,type="n",axes=F,xlab="",ylab="")
abline(v=0,lwd=lw,col="lightblue")
abline(h=yMean2,col="grey",lwd=lw,lty=2)
abline(h=yTen[1],col="grey",lwd=lw,lty=3)
abline(h=yTen[2],col="grey",lwd=lw,lty=3)
points(nucSelRealSCMTOjr$focal,vec2,col="red",cex=pointScaler)
par(new=T)
plot(nucSelRealSCMTOjr$focal,nucSelRealSCMTOjr$Nrate,ylim=yBounds1,type="n",axes=F,xlab="",ylab="")
axis(1,at=c(-125,125),labels=F,tick=T)
axis(1,at=c(-73,0,73))
points(nucSelRealSCMTOjr$focal,nucSelRealSCMTOjr$NormRate,col="darkblue",cex=pointScaler)
ss1<-smooth.spline(nucSelRealSCMTOjr$focal,nucSelRealSCMTOjr$NormRate,df=18)
points(ss1,lwd=lw,type="l",col="darkblue",cex=pointScaler)
axis(2,at=c(yBounds1[1],yMean1,yBounds1[2]),labels=round(c(yBounds1[1],yMean1,yBounds1[2]),digits=4),col="darkblue",col.axis="darkblue")
par(new=T)
plot(nucSelRealSCMTOjr$focal,vec2,ylim=yBounds2,type="n",axes=F,xlab="",ylab="")
ss<-smooth.spline(nucSelRealSCMTOjr$focal,vec2,df=18)
points(ss,lwd=3,type="l",col="red")
axis(4,at=c(yBounds2[1],yBounds2[2]),labels=c("     -15%","+15%"),tck=0.05,col="red",col.axis="red",mgp=c(0,-2,0))
axis(4,at=c(yBounds2[1],yMean2,yBounds2[2]),labels=round(c(yBounds2[1],yMean2,yBounds2[2]),digits=4),col="red",col.axis="red")
mtext("Okazaki junction rate (per nt)",side=2,line=3,col="darkblue")
mtext("Polymorphism rate (per nt)",side=4,line=3,col="red")
mtext("Distance from nucleosome dyad (nt)",side=1,line=3,col="black")

