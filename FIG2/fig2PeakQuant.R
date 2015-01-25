#############################
# 5 quantification
Cairo(file="fivepe.pdf",type="pdf",height=80,width=120,units="mm",bg="white",canvas="white",pointsize=10)

#backgroundWindow<-2001
#ojrWTrep$forDD2kbg<-rollapply(ojrWTrep$forCovDD,width=backgroundWindow,by=1,FUN=mean,align="center",na.rm=FALSE,fill=NA)
#ojrWTrep$revDD2kbg<-rollapply(ojrWTrep$revCovDD,width=backgroundWindow,by=1,FUN=mean,align="center",na.rm=FALSE,fill=NA)
#ojrWTrepAsymVec<-log2((ojrWTrep$forDD2kbg+1)/(ojrWTrep$revDD2kbg+1))


qthr<-quantile(ojrWTrep$forClobber[which((ojrWTrep$forDDbg[gm$gid]+ojrWTrep$revDDbg[gm$gid])>30 & gm$exclude[gm$gid]==0 & ojrWTrepAsymVec[gm$gid]>0)],probs=c(0.8,0.99,0.999),na.rm=T)


fr<-c(-5:5)
sigCall<-matrix(gm2$polCall,gidLn,1)
sigCount<-matrix(gm2$polCount,gidLn,1)
focalLeft<-c(1:5)
focalRight<-c(7:11)

ZP<-which((ojrWTrep$forDDbg[gm$gid]+ojrWTrep$revDDbg[gm$gid])>30 & gm$exclude[gm$gid]==0 & ojrWTrepAsymVec[gm$gid]>0 & ojrWTrep$forClobber[gm$gid] > as.numeric(qthr[3]))

nBoots<-100
lfts<-rep(0,nBoots)
rghts<-rep(0,nBoots)
for(i in 1:nBoots){
 ZPboot<-sample(ZP,size=length(ZP),replace=T)
 cmPolCount<-signalCounterMatrixTri(nSites=ZPboot,focalRange=fr,triTransform=gm2$gid,signal=sigCount,sqDF=gm2,mask=gm2$exclude)
 cmPolCall<-signalCounterMatrixTri(nSites=ZPboot,focalRange=fr,triTransform=gm2$gid,signal=sigCall,sqDF=gm2,mask=gm2$exclude)
 lfts[i]<-mean(cmPolCount$NormRate[focalLeft]/cmPolCall$NormRate[focalLeft])
 rghts[i]<-mean(cmPolCount$NormRate[focalRight]/cmPolCall$NormRate[focalRight])
 print(i)
}

cmPolCount<-signalCounterMatrixTri(nSites=ZPboot,focalRange=fr,triTransform=gm2$gid,signal=sigCount,sqDF=gm2,mask=gm2$exclude)
cmPolCall<-signalCounterMatrixTri(nSites=ZPboot,focalRange=fr,triTransform=gm2$gid,signal=sigCall,sqDF=gm2,mask=gm2$exclude)
lft<-mean(cmPolCount$NormRate[focalLeft]/cmPolCall$NormRate[focalLeft])
rght<-mean(cmPolCount$NormRate[focalRight]/cmPolCall$NormRate[focalRight])
lftKeep1<-lft
rghtKeep1<-rght
lftsKeep1<-lfts
rghtsKeep1<-rghts

for(i in 1:nBoots){
 iCol<-i+2
 cmPolCount<-signalCounterMatrixTri(nSites=ZP,focalRange=fr,triTransform=triShuff[,iCol],signal=sigCount,sqDF=gm2,mask=gm2$exclude)
 cmPolCall<-signalCounterMatrixTri(nSites=ZP,focalRange=fr,triTransform=triShuff[,iCol],signal=sigCall,sqDF=gm2,mask=gm2$exclude)
 lfts[i]<-mean(cmPolCount$NormRate[focalLeft]/cmPolCall$NormRate[focalLeft])
 rghts[i]<-mean(cmPolCount$NormRate[focalRight]/cmPolCall$NormRate[focalRight])
 print(paste("triQuant",i))
}
lftsTriKeep1<-lfts
rghtsTriKeep1<-rghts




ZP<-which((ojrWTrep$forDDbg[gm$gid]+ojrWTrep$revDDbg[gm$gid])>30 & gm$exclude[gm$gid]==0 & ojrWTrepAsymVec[gm$gid]>0 & ojrWTrep$forClobber[gm$gid] > as.numeric(qthr[2]) & ojrWTrep$forClobber[gm$gid] < as.numeric(qthr[3]))
xp<-3

nBoots<-100
lfts<-rep(0,nBoots)
rghts<-rep(0,nBoots)
for(i in 1:nBoots){
 ZPboot<-sample(ZP,size=length(ZP),replace=T)
 cmPolCount<-signalCounterMatrixTri(nSites=ZPboot,focalRange=fr,triTransform=gm2$gid,signal=sigCount,sqDF=gm2,mask=gm2$exclude)
 cmPolCall<-signalCounterMatrixTri(nSites=ZPboot,focalRange=fr,triTransform=gm2$gid,signal=sigCall,sqDF=gm2,mask=gm2$exclude)
 lfts[i]<-mean(cmPolCount$NormRate[focalLeft]/cmPolCall$NormRate[focalLeft])
 rghts[i]<-mean(cmPolCount$NormRate[focalRight]/cmPolCall$NormRate[focalRight])
 print(i)
}

cmPolCount<-signalCounterMatrixTri(nSites=ZPboot,focalRange=fr,triTransform=gm2$gid,signal=sigCount,sqDF=gm2,mask=gm2$exclude)
cmPolCall<-signalCounterMatrixTri(nSites=ZPboot,focalRange=fr,triTransform=gm2$gid,signal=sigCall,sqDF=gm2,mask=gm2$exclude)
lft<-mean(cmPolCount$NormRate[focalLeft]/cmPolCall$NormRate[focalLeft])
rght<-mean(cmPolCount$NormRate[focalRight]/cmPolCall$NormRate[focalRight])
lftKeep2<-lft
rghtKeep2<-rght
lftsKeep2<-lfts
rghtsKeep2<-rghts



for(i in 1:nBoots){
 iCol<-i+2
 cmPolCount<-signalCounterMatrixTri(nSites=ZP,focalRange=fr,triTransform=triShuff[,iCol],signal=sigCount,sqDF=gm2,mask=gm2$exclude)
 cmPolCall<-signalCounterMatrixTri(nSites=ZP,focalRange=fr,triTransform=triShuff[,iCol],signal=sigCall,sqDF=gm2,mask=gm2$exclude)
 lfts[i]<-mean(cmPolCount$NormRate[focalLeft]/cmPolCall$NormRate[focalLeft])
 rghts[i]<-mean(cmPolCount$NormRate[focalRight]/cmPolCall$NormRate[focalRight])
 print(paste("triQuant",i))
}
lftsTriKeep2<-lfts
rghtsTriKeep2<-rghts




validSites<-which((ojrWTrep$forDDbg[gm$gid]+ojrWTrep$revDDbg[gm$gid])>30 & gm$exclude[gm$gid]==0 & ojrWTrepAsymVec[gm$gid]>0 & !is.na(ojrWTrep$forClobber[gm$gid]) & ojrWTrep$forClobber[gm$gid] < as.numeric(qthr[2]))
ZP<-gm$gid[sample(validSites,size=100000,replace=T,prob=ojrWTrep$forClobber[validSites])]


nBoots<-100
lfts<-rep(0,nBoots)
rghts<-rep(0,nBoots)
for(i in 1:nBoots){
 ZPboot<-sample(ZP,size=100000,replace=T,prob=ojrWTrep$forClobber[ZP])
 cmPolCount<-signalCounterMatrixTri(nSites=ZPboot,focalRange=fr,triTransform=gm2$gid,signal=sigCount,sqDF=gm2,mask=gm2$exclude)
 cmPolCall<-signalCounterMatrixTri(nSites=ZPboot,focalRange=fr,triTransform=gm2$gid,signal=sigCall,sqDF=gm2,mask=gm2$exclude)
 lfts[i]<-mean(cmPolCount$NormRate[focalLeft]/cmPolCall$NormRate[focalLeft])
 rghts[i]<-mean(cmPolCount$NormRate[focalRight]/cmPolCall$NormRate[focalRight])
 print(i)
}

#cmPolCount<-signalCounterMatrixTri(nSites=ZPboot,focalRange=fr,triTransform=gm2$gid,signal=sigCount,sqDF=gm2,mask=gm2$exclude)
#cmPolCall<-signalCounterMatrixTri(nSites=ZPboot,focalRange=fr,triTransform=gm2$gid,signal=sigCall,sqDF=gm2,mask=gm2$exclude)
#lft<-mean(cmPolCount$NormRate[focalLeft]/cmPolCall$NormRate[focalLeft])
#rght<-mean(cmPolCount$NormRate[focalRight]/cmPolCall$NormRate[focalRight])
lft<-median(lfts)
lftKeep3<-lft
rghtKeep3<-median(rghts)
lftsKeep3<-lfts
rghtsKeep3<-rghts



for(i in 1:nBoots){
 iCol<-i+2
 cmPolCount<-signalCounterMatrixTri(nSites=ZP,focalRange=fr,triTransform=triShuff[,iCol],signal=sigCount,sqDF=gm2,mask=gm2$exclude)
 cmPolCall<-signalCounterMatrixTri(nSites=ZP,focalRange=fr,triTransform=triShuff[,iCol],signal=sigCall,sqDF=gm2,mask=gm2$exclude)
 lfts[i]<-mean(cmPolCount$NormRate[focalLeft]/cmPolCall$NormRate[focalLeft])
 rghts[i]<-mean(cmPolCount$NormRate[focalRight]/cmPolCall$NormRate[focalRight])
 print(paste("triQuant",i))
}
lftsTriKeep3<-lfts
rghtsTriKeep3<-rghts



xp<-1
###
# plotting it
qBounds<-c(0.025,0.975)
par(mar=c(4,4,2,2))
plot(c(0,9),c(0.02,0.04),type="n",axes=F,xlab="",ylab="")
#axis(1,at=c(1,2,3,4,5,6),labels=c("3'","5'","3'","5'","3'","5'"))
axis(1,at=c(0,1.5,4.5,7.5,9),labels=c("","Top 0.1%","Top 1%","prob. sampled",""))
axis(2,col="red",col.axis="red")
#mtext("Polymorphism rate",side=2,line=3,col="red")
# triShuff
#points(c(1,1),c(as.numeric(quantile(lftsTriKeep1,probs=qBounds))),type="l",col="pink")
meanPoint<-mean(lftsTriKeep1)
mpsd<-sd(lftsTriKeep1)
quant<-c(meanPoint-mpsd,meanPoint+mpsd)
#polygon(c(xp-.2,xp+.2,xp+.2,xp-.2),c(quant[1],quant[1],quant[2],quant[2]),border="pink")
points(c(xp+1,xp+1),c(quant[1],quant[2]),type="l",col="pink")
points(c(xp+1),c(meanPoint),col="pink",pch=20,cex=2)
medPoint<-median(lftsTriKeep1)
#points(c(0.8,1.2),c(medPoint,medPoint),col="pink",lwd=3,type="l",lend=1)
#points(c(2,2),c(as.numeric(quantile(rghtsTriKeep1,probs=qBounds))),type="l",col="pink")
#quant<-as.numeric(quantile(rghtsTriKeep1,probs=c(.33,.66)))
meanPoint<-mean(rghtsTriKeep1)
mpsd<-sd(rghtsTriKeep1)
quant<-c(meanPoint-mpsd,meanPoint+mpsd)
#polygon(c(xp+.8,xp+1.2,xp+1.2,xp+.8),c(quant[1],quant[1],quant[2],quant[2]),border="pink")
points(c(xp,xp),c(quant[1],quant[2]),type="l",col="pink")
points(c(xp),c(meanPoint),col="pink",pch=20,cex=2)
medPoint<-median(rghtsTriKeep1)
#points(c(1.8,2.2),c(medPoint,medPoint),col="pink",lwd=3,type="l",lend=1)
# real observations with bootstrap
lquant<-quantile(lftsKeep1,probs=c(.1,.9))
rquant<-quantile(rghtsKeep1,probs=c(.1,.9))
meanPoint<-mean(lftsKeep1)
mpsd<-sd(lftsKeep1)
quant<-c(meanPoint-mpsd,meanPoint+mpsd)
#polygon(c(xp-.2,xp+.2,xp+.2,xp-.2),c(quant[1],quant[1],quant[2],quant[2]),border="red")
points(c(xp,xp),c(quant[1],quant[2]),type="l",col="red")
points(c(xp),c(meanPoint),col="red",pch=20,cex=2)

#points(c(1,1),c(as.numeric(quantile(lftsKeep1,probs=qBounds))),type="l",col="red")
#points(c(0.8,1.2),c(lftKeep1,lftKeep1),col="red",lwd=3,type="l",lend=1)
meanPoint<-mean(rghtsKeep1)
mpsd<-sd(rghtsKeep1)
quant<-c(meanPoint-mpsd,meanPoint+mpsd)
#polygon(c(xp+.8,xp+1.2,xp+1.2,xp+.8),c(quant[1],quant[1],quant[2],quant[2]),border="red")
points(c(xp+1,xp+1),c(quant[1],quant[2]),type="l",col="red")
points(c(xp+1),c(meanPoint),col="red",pch=20,cex=2)
#points(c(2,2),c(as.numeric(quantile(rghtsKeep1,probs=qBounds))),type="l",col="red")
#points(c(1.8,2.2),c(rghtKeep1,rghtKeep1),col="red",lwd=3,type="l",lend=1)

xp<-4

#plot(c(0,6),c(0.02,0.03),type="n")
# triShuff
#points(c(xp,xp),c(as.numeric(quantile(lftsTriKeep2,probs=qBounds))),type="l",col="pink")
meanPoint<-mean(lftsTriKeep2)
mpsd<-sd(lftsTriKeep2)
quant<-c(meanPoint-mpsd,meanPoint+mpsd)
#polygon(c(xp-.2,xp+.2,xp+.2,xp-.2),c(quant[1],quant[1],quant[2],quant[2]),border="pink")
lcol<-"pink"
lmod<-0
points(c(xp+lmod,xp+lmod),c(quant[1],quant[2]),type="l",col=lcol)
points(c(xp+lmod),c(meanPoint),col=lcol,pch=20,cex=2)

medPoint<-median(lftsTriKeep2)
#points(c(xp-.2,xp+.2),c(medPoint,medPoint),col="pink",lwd=3,type="l",lend=1)
#points(c(xp+1,xp+1),c(as.numeric(quantile(rghtsTriKeep2,probs=qBounds))),type="l",col="pink")
meanPoint<-mean(rghtsTriKeep2)
mpsd<-sd(rghtsTriKeep2)
quant<-c(meanPoint-mpsd,meanPoint+mpsd)
#polygon(c(xp+.8,xp+1.2,xp+1.2,xp+.8),c(quant[1],quant[1],quant[2],quant[2]),border="pink")
lcol<-"pink"
lmod<-1
points(c(xp+lmod,xp+lmod),c(quant[1],quant[2]),type="l",col=lcol)
points(c(xp+lmod),c(meanPoint),col=lcol,pch=20,cex=2)

medPoint<-median(rghtsTriKeep2)
#points(c(xp+.8,xp+1.2),c(medPoint,medPoint),col="pink",lwd=3,type="l",lend=1)
# real observations with bootstrap

lquant<-quantile(lftsKeep2,probs=c(.1,.9))
rquant<-quantile(rghtsKeep2,probs=c(.1,.9))
meanPoint<-mean(lftsKeep2)
mpsd<-sd(lftsKeep2)
quant<-c(meanPoint-mpsd,meanPoint+mpsd)
#polygon(c(xp-.2,xp+.2,xp+.2,xp-.2),c(quant[1],quant[1],quant[2],quant[2]),border="red")
lcol<-"red"
lmod<-0
points(c(xp+lmod,xp+lmod),c(quant[1],quant[2]),type="l",col=lcol)
points(c(xp+lmod),c(meanPoint),col=lcol,pch=20,cex=2)

#points(c(xp,xp),c(as.numeric(quantile(lftsKeep2,probs=qBounds))),type="l",col="red")
#points(c(xp-.2,xp+.2),c(lftKeep2,lftKeep2),col="red",lwd=3,type="l",lend=1)
meanPoint<-mean(rghtsKeep2)
mpsd<-sd(rghtsKeep2)
quant<-c(meanPoint-mpsd,meanPoint+mpsd)
#polygon(c(xp+.8,xp+1.2,xp+1.2,xp+.8),c(quant[1],quant[1],quant[2],quant[2]),border="red")
lcol<-"red"
lmod<-1
points(c(xp+lmod,xp+lmod),c(quant[1],quant[2]),type="l",col=lcol)
points(c(xp+lmod),c(meanPoint),col=lcol,pch=20,cex=2)

#points(c(xp+1,xp+1),c(as.numeric(quantile(rghtsKeep2,probs=qBounds))),type="l",col="red")
#points(c(xp+.8,xp+1.2),c(rghtKeep2,rghtKeep2),col="red",lwd=3,type="l",lend=1)

xp<-7

#plot(c(0,6),c(0.02,0.03),type="n")
# triShuff
#points(c(xp,xp),c(as.numeric(quantile(lftsTriKeep3,probs=qBounds))),type="l",col="pink")
meanPoint<-mean(lftsTriKeep3)
mpsd<-sd(lftsTriKeep3)
quant<-c(meanPoint-mpsd,meanPoint+mpsd)
#polygon(c(xp-.2,xp+.2,xp+.2,xp-.2),c(quant[1],quant[1],quant[2],quant[2]),border="pink")
lcol<-"pink"
lmod<-0
points(c(xp+lmod,xp+lmod),c(quant[1],quant[2]),type="l",col=lcol)
points(c(xp+lmod),c(meanPoint),col=lcol,pch=20,cex=2)

medPoint<-median(lftsTriKeep3)
#points(c(xp-.2,xp+.2),c(medPoint,medPoint),col="pink",lwd=3,type="l",lend=1)
#points(c(xp+1,xp+1),c(as.numeric(quantile(rghtsTriKeep3,probs=qBounds))),type="l",col="pink")
meanPoint<-mean(rghtsTriKeep3)
mpsd<-sd(rghtsTriKeep3)
quant<-c(meanPoint-mpsd,meanPoint+mpsd)
#polygon(c(xp+.8,xp+1.2,xp+1.2,xp+.8),c(quant[1],quant[1],quant[2],quant[2]),border="pink")
lcol<-"pink"
lmod<-1
points(c(xp+lmod,xp+lmod),c(quant[1],quant[2]),type="l",col=lcol)
points(c(xp+lmod),c(meanPoint),col=lcol,pch=20,cex=2)

medPoint<-median(rghtsTriKeep3)
#points(c(xp+.8,xp+1.2),c(medPoint,medPoint),col="pink",lwd=3,type="l",lend=1)
# real observations with bootstrap

lquant<-quantile(lftsKeep3,probs=c(.1,.9))
rquant<-quantile(rghtsKeep3,probs=c(.1,.9))
meanPoint<-mean(lftsKeep3)
mpsd<-sd(lftsKeep3)
quant<-c(meanPoint-mpsd,meanPoint+mpsd)
#polygon(c(xp-.2,xp+.2,xp+.2,xp-.2),c(quant[1],quant[1],quant[2],quant[2]),border="red")
lcol<-"red"
lmod<-0
points(c(xp+lmod,xp+lmod),c(quant[1],quant[2]),type="l",col=lcol)
points(c(xp+lmod),c(meanPoint),col=lcol,pch=20,cex=2)

#points(c(xp,xp),c(as.numeric(quantile(lftsKeep3,probs=qBounds))),type="l",col="red")
#points(c(xp-.2,xp+.2),c(lftKeep3,lftKeep3),col="red",lwd=3,type="l",lend=1)
meanPoint<-mean(rghtsKeep3)
mpsd<-sd(rghtsKeep3)
quant<-c(meanPoint-mpsd,meanPoint+mpsd)
#polygon(c(xp+.8,xp+1.2,xp+1.2,xp+.8),c(quant[1],quant[1],quant[2],quant[2]),border="red")
lcol<-"red"
lmod<-1
points(c(xp+lmod,xp+lmod),c(quant[1],quant[2]),type="l",col=lcol)
points(c(xp+lmod),c(meanPoint),col=lcol,pch=20,cex=2)

#points(c(xp+1,xp+1),c(as.numeric(quantile(rghtsKeep3,probs=qBounds))),type="l",col="red")
#points(c(xp+.8,xp+1.2),c(rghtKeep3,rghtKeep3),col="red",lwd=3,type="l",lend=1)

dev.off()


