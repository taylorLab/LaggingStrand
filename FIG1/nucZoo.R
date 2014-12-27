###
# nucleosome zoo rate
fr<-c(-125:125)
sig<-matrix(gm$zoo,gidLn,1)
nucSelRealSCMTzoo<-signalCounterMatrixTri(nSites=nucSel,focalRange=fr,triTransform=gm$gid,signal=sig,sqDF=gm,mask=gm$exclude)
#sig<-matrix(c((ojrWTrep$forClobber+ojrWTrep$revClobber)/2,gidLn,1))
#nucSelRealSCMTOjr<-signalCounterMatrixTri(nSites=nucSel,focalRange=fr,triTransform=gm$gid,signal=sig,sqDF=gm,mask=gm$exclude)

sig<-matrix(gm$polCall,gidLn,1)
nucSelRealSCMTpolCall<-signalCounterMatrixTri(nSites=nucSel,focalRange=fr,triTransform=gm$gid,signal=sig,sqDF=gm,mask=gm$exclude)
sig<-matrix(gm$polCount,gidLn,1)
nucSelRealSCMTpol<-signalCounterMatrixTri(nSites=nucSel,focalRange=fr,triTransform=gm$gid,signal=sig,sqDF=gm,mask=gm$exclude)

vec1<-c(nucSelRealSCMTpol$NormRate/nucSelRealSCMTpolCall$NormRate)
vec2<-nucSelRealSCMTzoo$Nrate


ss2<-smooth.spline(fr,vec2,df=18)
ss2p<-predict(ss2,fr)
ss1<-smooth.spline(fr,vec1,df=18)
ss1p<-predict(ss1,fr)

yMean1<-mean(ss1p$y)
yMean2<-mean(ss2p$y)
yTen<-c(yMean2*.98,yMean2*1.02)
yTen1<-c(yMean1*.9,yMean1*1.1)

yMarg=.15
yBounds1<-c(yMean1*(1-yMarg),yMean1*(1+yMarg))
yMarg2=0.03
yBounds2<-c(yMean2*(1-yMarg2),yMean2*(1+yMarg2))



plot(ss2,ylim=yBounds2,type="n",axes=F,xlab="",ylab="")
abline(v=0,lwd=lw,col="lightblue")
abline(h=yMean2,col="grey",lwd=lw,lty=2)
#abline(h=yTen[1],col="grey",lwd=3,lty=3)
#abline(h=yTen[2],col="grey",lwd=3,lty=3)
#points(fr,vec2,col="orange",cex=0.5)
points(ss2,lwd=lw,type="l",col=lineSubr)
par(new=T)
plot(ss1,ylim=yBounds1,type="n",axes=F,xlab="",ylab="")
axis(1,at=c(-125,125),labels=F,tick=T)
axis(1,at=c(-73,0,73))
points(ss1,lwd=lw,type="l",col="red")
axis(2,at=c(yTen1[1],yMean1,yTen1[2]),labels=c("-10%","mean","+10%"),col="red",col.axis="red")
par(new=T)
plot(ss2,ylim=yBounds2,type="n",axes=F,xlab="",ylab="")
axis(4,at=c(yTen[1],yMean2,yTen[2]),labels=c("-2%","mean","+2%"),col=lineSubr,col.axis=lineSubr)
#mtext("Relative sub. rate",side=4,line=3,col=lineSubr)
mtext("Nucleosome (nt)",side=1,line=3,col="black")
