
###
# reb zoo rate
sig<-matrix(gm$zoo,gidLn,1)
fr<-c(-80:80)
rebSelRealSCMTzoo<-signalCounterMatrixTri(nSites=rebSel,focalRange=fr,triTransform=gm$gid,signal=sig,sqDF=gm,mask=gm$exclude)
sig<-matrix(c((ojrWTrep$forClobber+ojrWTrep$revClobber)/2,gidLn,1))
rebSelRealSCMTOjr<-signalCounterMatrixTri(nSites=rebSel,focalRange=fr,triTransform=gm$gid,signal=sig,sqDF=gm,mask=gm$exclude)

sig<-matrix(gm$zoo,gidLn,1)
rebSelRealSCMTzoo<-signalCounterMatrixTri(nSites=rebSel,focalRange=fr,triTransform=gm$gid,signal=sig,sqDF=gm,mask=gm$exclude)

sig<-matrix(gm2$polCall,gidLn,1)
rebSelRealSCMTpolCall<-signalCounterMatrixTri(nSites=rebSel,focalRange=fr,triTransform=gm$gid,signal=sig,sqDF=gm,mask=gm$exclude)
sig<-matrix(gm2$polCount,gidLn,1)
rebSelRealSCMTpol<-signalCounterMatrixTri(nSites=rebSel,focalRange=fr,triTransform=gm$gid,signal=sig,sqDF=gm,mask=gm$exclude)


vec1<-c(rebSelRealSCMTpol$NormRate/rebSelRealSCMTpolCall$NormRate)
vec2<-c(rebSelRealSCMTzoo$NormRate)

ss2<-smooth.spline(fr,vec2,df=38)
ss2p<-predict(ss2,fr)
ss1<-smooth.spline(fr,vec1,df=38)
ss1p<-predict(ss1,fr)


frl<-length(fr)
frlLeft<-(frl/2)-20
frlRight<-(frl/2)+20
meanVec<-round(c(1:frlLeft,frlRight:frl))
yMean1<-mean(ss1p$y[meanVec])
yMean2<-mean(ss2p$y[meanVec])
yTen<-c(yMean2*.8,yMean2*1.2)
yTen1<-c(yMean1*.5,yMean1*1.5)


yMarg=1.2
yBounds1<-c(yMean1*(1-yMarg),yMean1*(1+yMarg))
yMarg2=0.35
yBounds2<-c(yMean2*(1-yMarg2),yMean2*(1+yMarg2))




plot(ss2,ylim=yBounds2,type="n",axes=F,xlab="",ylab="")
abline(v=0,lwd=lw,col="lightblue")
abline(h=yMean2,col="grey",lwd=lw,lty=2)
#abline(h=yTen[1],col="grey",lwd=3,lty=3)
#abline(h=yTen[2],col="grey",lwd=3,lty=3)
#points(rebSelRealSCMTOjr$focal,vec2,col="orange",cex=0.5)
points(ss2,lwd=lw,type="l",col=lineSubr)

par(new=T)
plot(ss1,ylim=yBounds1,type="n",axes=F,xlab="",ylab="")
axis(1,at=c(min(fr),max(fr)),labels=F,tick=T)
axis(1,at=c(-50,0,50))
points(ss1,lwd=lw,type="l",col="red")
#axis(2,at=round(c(yBounds1[1],yMean1,yBounds1[2]),digits=4),col="darkblue",col.axis="darkblue")
axis(2,at=c(yMean1*.2,yMean1,yMean1*2),labels=c("-80%","mean","+100%"),col="red",col.axis="red")

par(new=T)
plot(ss2,ylim=yBounds2,type="n",axes=F,xlab="",ylab="")
#points(ss2,lwd=4,type="l",col="orange")
#axis(4,at=round(c(yBounds2[1],yMean2,yBounds2[2]),digits=4),col="orange",col.axis="orange")
axis(4,at=round(c(yTen[1],yMean2,yTen[2]),digits=4),labels=c("-20%","mean","+20%"),col=lineSubr,col.axis=lineSubr)

#mtext("Okazaki junction rate (per nt)",side=2,line=3,col="darkblue")
mtext("Relative sub. rate",side=4,line=3,col=lineSubr)
mtext("Reb1 site (nt)",side=1,line=3,col="black")

