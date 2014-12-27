
rebSel<-reb$gid[which((reb$class=="primary" | reb$class=="monomer") & (ojrWTrep$forDDbg[reb$gid]+ojrWTrep$revDDbg[reb$gid])>30)]

fr<-c(-80:80)
sig<-matrix(c((ojrWTrep$forClobber+ojrWTrep$revClobber)/2,gidLn,1))
rebSelRealSCMTOjr<-signalCounterMatrixTri(nSites=rebSel,focalRange=fr,triTransform=gm$gid,signal=sig,sqDF=gm,mask=gm$exclude)
sig<-matrix(gm2$polCall,gidLn,1)
rebSelRealSCMTpolCall<-signalCounterMatrixTri(nSites=rebSel,focalRange=fr,triTransform=gm$gid,signal=sig,sqDF=gm,mask=gm$exclude)
sig<-matrix(gm2$polCount,gidLn,1)
rebSelRealSCMTpol<-signalCounterMatrixTri(nSites=rebSel,focalRange=fr,triTransform=gm$gid,signal=sig,sqDF=gm,mask=gm$exclude)


vec2<-c(rebSelRealSCMTpol$NormRate/rebSelRealSCMTpolCall$NormRate)
frl<-length(fr)
frlLeft<-(frl/2)-20
frlRight<-(frl/2)+20
meanVec<-round(c(1:frlLeft,frlRight:frl))
yMean1<-mean(rebSelRealSCMTOjr$NormRate[meanVec])
yMean2<-mean(vec2[meanVec])
yTen<-c(yMean2*.9,yMean2*1.1)
yMarg<-c(max(rebSelRealSCMTOjr$NormRate/yMean1,yMean1/rebSelRealSCMTOjr$NormRate,vec2/yMean2,yMean2/vec2))-1
yBounds1<-c(yMean1*.2,yMean1*2.5)
yBounds2<-c(yMean2*.2,yMean2*2.5)

plot(rebSelRealSCMTOjr$focal,vec2,ylim=yBounds2,type="n",axes=F,xlab="",ylab="")
abline(v=0,lwd=lw,col="lightblue")
abline(h=yMean2,col="grey",lwd=lw,lty=2)
abline(h=yTen[1],col="grey",lwd=lw,lty=3)
abline(h=yTen[2],col="grey",lwd=lw,lty=3)
points(rebSelRealSCMTOjr$focal,vec2,col="red",cex=pointScaler)
par(new=T)
plot(rebSelRealSCMTOjr$focal,rebSelRealSCMTOjr$NormRate,ylim=yBounds1,type="n",axes=F,xlab="",ylab="")
axis(1,at=c(min(fr),max(fr)),labels=F,tick=T)
axis(1,at=c(-50,0,50))
ss1<-smooth.spline(rebSelRealSCMTOjr$focal,rebSelRealSCMTOjr$Nrate,df=48)
points(ss1,lwd=3,type="l",col="darkblue")

points(rebSelRealSCMTOjr$focal,rebSelRealSCMTOjr$NormRate,col="darkblue",cex=pointScaler)
axis(2,at=round(c(yBounds1[1],yMean1,yBounds1[2]),digits=4),col="darkblue",col.axis="darkblue")
par(new=T)
plot(rebSelRealSCMTOjr$focal,vec2,ylim=yBounds2,type="n",axes=F,xlab="",ylab="")
ss<-smooth.spline(rebSelRealSCMTOjr$focal,vec2,df=48)
points(ss,lwd=3,type="l",col="red")
axis(4,at=c(yBounds2[1],yMean2*2),labels=c("     -80%","+100%"),tck=0.05,col="red",col.axis="red",mgp=c(0,-2.5,0))
axis(4,at=round(c(yBounds2[1],yMean2,yBounds2[2]),digits=4),col="red",col.axis="red")
mtext("Okazaki junction rate (per nt)",side=2,line=3,col="darkblue")
mtext("Polymorphism rate (per nt)",side=4,line=3,col="red")
mtext("Distance from Reb1 binding site (nt)",side=1,line=3,col="black")
