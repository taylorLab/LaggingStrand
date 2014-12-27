####
# triComp null for nucleosomes
dateStamp<-format(Sys.time(),"%Y%m%d")
fr<-c(-125:125)
dfree<-18
nucf30o80c30<-nucs$gid[which(nucs$occupancy>80 & nucs$fuzzyness<30 & nucs$tssProximity==0 & (ojrWTrep$forCovDD[nucs$gid]+ojrWTrep$revCovDD[nucs$gid])>30)]
nucSel<-nucf30o80c30
nTri<-100
if(!exists("nmat1")){
 nmat1<-matrix(rep(0,length(fr)*nTri),nTri,length(fr))
 nmat2<-nmat1
 for(i in 1:nTri){
  iCol<-i+2
  sig<-matrix(c((ojrWTrep$forClobber+ojrWTrep$revClobber)/2,gidLn,1))
  nucSelRealSCMTOjrTri1<-signalCounterMatrixTri(nSites=nucSel,focalRange=fr,triTransform=triShuff[,iCol],signal=sig,sqDF=gm,mask=gm$exclude)
  sig<-matrix(gm2$polCall,gidLn,1)
  nucSelRealSCMTpolCallTri1<-signalCounterMatrixTri(nSites=nucSel,focalRange=fr,triTransform=triShuff[,iCol],signal=sig,sqDF=gm,mask=gm$exclude)
  sig<-matrix(gm2$polCount,gidLn,1)
  nucSelRealSCMTpolTri1<-signalCounterMatrixTri(nSites=nucSel,focalRange=fr,triTransform=triShuff[,iCol],signal=sig,sqDF=gm,mask=gm$exclude)
  nmat1[i,]=predict(smooth.spline(fr,nucSelRealSCMTOjrTri1$NormRate,df=dfree),fr)$y
  nmat2[i,]=predict(smooth.spline(fr,c(nucSelRealSCMTpolTri1$NormRate/nucSelRealSCMTpolCallTri1$NormRate),df=dfree))$y
 }
}
write.table(nmat1,file=paste("nucTri.nmat1.",dateStamp,".Rtable"),quote=F,sep='\t',row.names=F,col.names=F)
write.table(nmat2,file=paste("nucTri.nmat2.",dateStamp,".Rtable"),quote=F,sep='\t',row.names=F,col.names=F)
mat1Margins<-apply(nmat1,2,quantile,probs=c(0.025,0.975))
mat1Means<-colMeans(nmat1)
meanVec<-round(c(1:frlLeft,frlRight:frl))
yMean1<-mean(mat1Means[meanVec])
yBounds1<-c(yMean1*.85,yMean1*1.15)
mat2Margins<-apply(nmat2,2,quantile,probs=c(0.025,0.975))
mat2Means<-colMeans(nmat2)
mat2TopMarginSSy<-mat2Margins[2,]
mat2BotMarginSSy<-mat2Margins[1,]

meanVec<-round(c(1:frlLeft,frlRight:frl))
yMean2<-mean(mat2Means[meanVec])
yBounds2<-c(yMean2*.85,yMean2*1.15)
yTen2<-c(yMean2*.9,yMean2*1.1)



plot(fr,nmat2[1,],ylim=yBounds2,type="n",axes=F,xlab="",ylab="")
points(fr,mat2Means,type="l",col="darkblue")
polygon(c(fr,rev(fr)),c(mat2TopMarginSSy,rev(mat2BotMarginSSy)),col="pink",border=F)

abline(v=0,lwd=lw,col="lightblue")
abline(h=yMean2,col="grey",lwd=lw,lty=2)
abline(h=yTen2[1],col="grey",lwd=lw,lty=3)
abline(h=yTen2[2],col="grey",lwd=lw,lty=3)
#axis(2,at=round(c(yBounds2[1],yMean2,yBounds2[2]),digits=3),col="red",col.axis="red")
axis(2,at=c(yBounds2[1],yMean2,yBounds2[2]),labels=c("-15%","mean","+15%"),col="red",col.axis="red")

ss2<-smooth.spline(fr,mat2Means,df=dfree)
points(ss2,type="l",col="red",lwd=lw)
par(new=T)
plot(fr,nmat1[1,],ylim=yBounds1,type="n",axes=F,xlab="",ylab="")
ss1<-smooth.spline(fr,mat1Means,df=dfree)
points(ss1,lwd=lw,type="l",col="darkblue")

mtext("Pol. rate",side=2,line=3,col="red")
mtext("Nucleosome (nt)",side=1,line=3,col="black")
axis(1,at=c(-125,125),labels=F,tick=T)
axis(1,at=c(-73,0,73))
