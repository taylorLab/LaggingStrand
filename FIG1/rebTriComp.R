###
# Reb1 compositional match

rebSel<-reb$gid[which((reb$class=="primary" | reb$class=="monomer") & (ojrWTrep$forDDbg[reb$gid]+ojrWTrep$revDDbg[reb$gid])>30)]
dateStamp<-format(Sys.time(),"%Y%m%d")

fr<-c(-80:80)
frl<-length(fr)
frlLeft<-(frl/2)-20
frlRight<-(frl/2)+20
nTri<-100
if(!exists("mat1")){
 mat1<-matrix(rep(0,length(fr)*nTri),nTri,length(fr))
 mat2<-mat1
 for(i in 1:nTri){
  iCol<-i+2
  sig<-matrix(c((ojrWTrep$forClobber+ojrWTrep$revClobber)/2,gidLn,1))
  rebSelRealSCMTOjrTri1<-signalCounterMatrixTri(nSites=rebSel,focalRange=fr,triTransform=triShuff[,iCol],signal=sig,sqDF=gm,mask=gm$exclude)
  sig<-matrix(gm2$polCall,gidLn,1)
  rebSelRealSCMTpolCallTri1<-signalCounterMatrixTri(nSites=rebSel,focalRange=fr,triTransform=triShuff[,iCol],signal=sig,sqDF=gm,mask=gm$exclude)
  sig<-matrix(gm2$polCount,gidLn,1)
  rebSelRealSCMTpolTri1<-signalCounterMatrixTri(nSites=rebSel,focalRange=fr,triTransform=triShuff[,iCol],signal=sig,sqDF=gm,mask=gm$exclude)

  mat1[i,]=predict(smooth.spline(fr,rebSelRealSCMTOjrTri1$NormRate,df=48),fr)$y
  mat2[i,]=predict(smooth.spline(fr,c(rebSelRealSCMTpolTri1$NormRate/rebSelRealSCMTpolCallTri1$NormRate),df=48))$y
 }
}
write.table(mat1,file=paste("rebTri.mat1.",dateStamp,".Rtable",sep=""),quote=F,sep='\t',row.names=F,col.names=F)
write.table(mat2,file=paste("rebTri.mat2.",dateStamp,".Rtable",sep=""),quote=F,sep='\t',row.names=F,col.names=F)
mat1Margins<-apply(mat1,2,quantile,probs=c(0.025,0.975))
mat1Means<-colMeans(mat1)
meanVec<-round(c(1:frlLeft,frlRight:frl))
yMean1<-mean(mat1Means[meanVec])
yBounds1<-c(yMean1*.2,yMean1*2.5)
mat2Margins<-apply(mat2,2,quantile,probs=c(0.1,0.9))
mat2TopMarginSSy<-mat2Margins[2,]
mat2BotMarginSSy<-mat2Margins[1,]
mat2Means<-colMeans(mat2)
meanVec<-round(c(1:frlLeft,frlRight:frl))
yMean2<-mean(mat2Means[meanVec])
yBounds2<-c(yMean2*.2,yMean2*2.5)
yTen2<-c(yMean2*.9,yMean2*1.1)

plot(fr,mat2Means,ylim=yBounds2,type="n",axes=F,xlab="",ylab="")
points(fr,mat2Means,type="l",col="darkblue")
polygon(c(fr,rev(fr)),c(mat2TopMarginSSy,rev(mat2BotMarginSSy)),col="pink",border=F)
abline(v=0,lwd=lw,col="lightblue")
abline(h=yMean2,col="grey",lwd=lw,lty=2)
#abline(h=yTen2[1],col="grey",lwd=3,lty=3)
#abline(h=yTen2[2],col="grey",lwd=3,lty=3)
#points(fr,mat2Means,type="l",col="red",lwd=3)
ss2<-smooth.spline(fr,mat2Means,df=48)
points(ss2,type="l",col="red",lwd=1)
axis(2,at=c(yBounds2[1],yMean2,yMean2*2),labels=c("-80%","mean","+100%"),col="red",col.axis="red")
par(new=T)
plot(fr,mat2Means,ylim=yBounds1,type="n",axes=F,xlab="",ylab="")
#points(fr,mat1Means,type="l",col="darkblue")
ss1<-smooth.spline(fr,mat1Means,df=48)
points(ss1,lwd=lw,type="l",col="darkblue")
#polygon(c(fr,rev(fr)),c(mat1Margins[1,],rev(mat1Margins[2,])),col="lightblue",border=F)
#mtext("Pol. rate",side=2,line=3,col="red")
mtext("Reb1 site (nt)",side=1,line=3,col="black")
axis(1,at=c(min(fr),max(fr)),labels=F,tick=T)
axis(1,at=c(-50,0,50))
