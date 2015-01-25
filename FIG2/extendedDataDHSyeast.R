
Cairo(file="suppFigDNase.pdf",type="pdf",height=80,width=220,units="mm",bg="white",canvas="white",pointsize=10)


layout(matrix(c(1,2,1,2,3,3,3,3,4,4,4,4),2,6))
#layout.show(4)

lw=1

####
# reb both strands

fr<-c(-80:80)

rebSel<-reb$gid[which((reb$class=="primary" | reb$class=="monomer") & (ojrWTrep$forDDbg[reb$gid]+ojrWTrep$revDDbg[reb$gid])>30 & gm$exclude[reb$gid]==0)]

ZP<-rebSel

sig<-matrix(footprintCoverage,gidLn,1)
cmFootprint<-signalCounterMatrixTri(nSites=ZP,focalRange=fr,triTransform=gm$gid,signal=sig,sqDF=gm,mask=gm$exclude)

sig<-matrix(ojrWTrep$forClobber,gidLn,1)
cmOjrFor<-signalCounterMatrixTri(nSites=ZP,focalRange=fr,triTransform=gm$gid,signal=sig,sqDF=gm,mask=gm$exclude)

sig<-matrix(ojrWTrep$revClobber,gidLn,1)
cmOjrRev<-signalCounterMatrixTri(nSites=ZP,focalRange=fr,triTransform=gm$gid,signal=sig,sqDF=gm,mask=gm$exclude)

sig<-matrix(ojrWTsmp$forClobber,gidLn,1)
cmOjrForSmp<-signalCounterMatrixTri(nSites=ZP,focalRange=fr,triTransform=gm$gid,signal=sig,sqDF=gm,mask=gm$exclude)

sig<-matrix(gm$zoo,gidLn,1)
cmZoo<-signalCounterMatrixTri(nSites=ZP,focalRange=fr,triTransform=gm$gid,signal=sig,sqDF=gm,mask=gm$exclude)

sig<-matrix(gm2$polCount,gidLn,1)
cmPolCount<-signalCounterMatrixTri(nSites=ZP,focalRange=fr,triTransform=gm$gid,signal=sig,sqDF=gm,mask=gm$exclude)

sig<-matrix(gm2$polCal,gidLn,1)
cmPolCall<-signalCounterMatrixTri(nSites=ZP,focalRange=fr,triTransform=gm$gid,signal=sig,sqDF=gm,mask=gm$exclude)


layout(matrix(c(1,2,1,2,3,3,3,3,4,4,4,4),2,6))

## reb two ways
par(mar=c(0,5,4,3))
yBounds<-c(0,max(cmFootprint$Nsignal))
plot(fr,cmFootprint$Nsignal,type="n",axes=F,xlab="",ylab="",ylim=yBounds)
polygon(c(fr,max(fr),min(fr)),c(cmFootprint$Nsignal,0,0),col="grey",border=NA)
abline(v=0,col="lightblue")
par(new=T)
plot(fr,cmOjrFor$NormRate,col="darkblue",type="l",lwd=lw,axes=F,xlab="",ylab="",ylim=c(0.004,0.015))
#arrows(-70,0.012,-40,0.012)
axis(2,col="darkblue",at=c(0.006,0.014),col.axis="darkblue")
par(mar=c(4,5,0,3))
plot(fr,cmFootprint$Nsignal,type="n",axes=F,xlab="",ylab="",ylim=rev(yBounds))
polygon(c(fr,max(fr),min(fr)),c(cmFootprint$Nsignal,0,0),col="grey",border=NA)
abline(v=0,col="lightblue")
par(new=T)
plot(fr,cmOjrRev$NormRate,col="darkblue",lwd=lw,type="l",axes=F,xlab="",ylab="",ylim=c(0.015,0.004))
#arrows(70,0.012,40,0.012)
axis(2,col="darkblue",at=c(0.006,0.014),col.axis="darkblue")
axis(1)
mtext("Reb1 site",side=1,line=2.5,col="black")
mtext("OJ rate (per nt)",side=2,line=3,col="darkblue")

## footprint mid


candMid<-dnase$V5+round((dnase$V6-dnase$V5)/2)
ZPmid<-candMid[which(gm$exclude[candMid]==0)]
ZP<-dnase$V5[which(gm$exclude[dnase$V5]==0)]

ZP<-ZPmid

sig<-matrix(gm$zoo,gidLn,1)
cmZoo<-signalCounterMatrixTri(nSites=ZP,focalRange=fr,triTransform=gm$gid,signal=sig,sqDF=gm,mask=gm$exclude)
cmZoor<-signalCounterMatrixTri(nSites=ZPr,focalRange=fr,triTransform=gm$gid,signal=sig,sqDF=gm,mask=gm$exclude)


sig<-matrix(footprintCoverage,gidLn,1)
cmFootprint<-signalCounterMatrixTri(nSites=ZP,focalRange=fr,triTransform=gm$gid,signal=sig,sqDF=gm,mask=gm$exclude)
cmFootprintr<-signalCounterMatrixTri(nSites=ZPr,focalRange=fr,triTransform=gm$gid,signal=sig,sqDF=gm,mask=gm$exclude)

sig<-matrix(ojrWTrep$forClobber,gidLn,1)
cmOjrFor<-signalCounterMatrixTri(nSites=ZP,focalRange=fr,triTransform=gm$gid,signal=sig,sqDF=gm,mask=gm$exclude)
cmOjrForr<-signalCounterMatrixTri(nSites=ZPr,focalRange=fr,triTransform=gm$gid,signal=sig,sqDF=gm,mask=gm$exclude)

sig<-matrix(ojrWTrep$revClobber,gidLn,1)
cmOjrRev<-signalCounterMatrixTri(nSites=ZP,focalRange=fr,triTransform=gm$gid,signal=sig,sqDF=gm,mask=gm$exclude)
cmOjrRevr<-signalCounterMatrixTri(nSites=ZPr,focalRange=fr,triTransform=gm$gid,signal=sig,sqDF=gm,mask=gm$exclude)


sig<-matrix(gm2$polCount,gidLn,1)
cmPolCount<-signalCounterMatrixTri(nSites=ZP,focalRange=fr,triTransform=gm$gid,signal=sig,sqDF=gm,mask=gm$exclude)
cmPolCountr<-signalCounterMatrixTri(nSites=ZPr,focalRange=fr,triTransform=gm$gid,signal=sig,sqDF=gm,mask=gm$exclude)

sig<-matrix(gm2$polCal,gidLn,1)
cmPolCall<-signalCounterMatrixTri(nSites=ZP,focalRange=fr,triTransform=gm$gid,signal=sig,sqDF=gm,mask=gm$exclude)
cmPolCallr<-signalCounterMatrixTri(nSites=ZPr,focalRange=fr,triTransform=gm$gid,signal=sig,sqDF=gm,mask=gm$exclude)

par(mar=c(4,5,4,5))
plot(fr,cmFootprint$Nsignal,type="n",axes=F,xlab="",ylab="")
axis(1)
polygon(c(fr,max(fr),min(fr)),c((cmFootprint$Nsignal),0,0),col="grey",border=NA)
#abline(v=0,col="lightblue")
par(new=T)
#plot(fr,cmOjrFor$NormRate,col="darkblue",type="l",lwd=lw,axes=F,xlab="",ylab="")
plot(fr,(cmOjrFor$NormRate+cmOjrRev$NormRate)/2,col="darkblue",type="l",lwd=lw,axes=F,xlab="",ylab="")
axis(2,col="darkblue",col.axis="darkblue")
par(new=T)
ssPol<-smooth.spline(fr,cmPolCount$NormRate/cmPolCall$NormRate,df=35)
plot(ssPol,col="red",lwd=lw,type="l",axes=F,xlab="",ylab="")
axis(4,col="red",col.axis="red")
mtext("DNase footprint mid (nt)",side=1,line=2.5,col="black")
mtext("Okazaki junction rate (per nt)",side=2,line=3,col="darkblue")
mtext("Polymorphism rate",side=4,line=3,col="red")


####
# ZPmid with Reb1 and Rap1 masking
rebRapMarginMask<-rep(0,gidLn)
for(i in 1:length(reb$gid)){
 rebRapMarginMask[(reb$gid[i]-50):(reb$gid[i]+50)]=1
}
for(i in 1:length(rap$V4)){
 rebRapMarginMask[(rap$V5[i]-50):(rap$V5[i]+50)]=1
}
ZPmidNoRapReb<-ZPmid[which(rebRapMarginMask[ZPmid]==0)]
ZP<-ZPmidNoRapReb


sig<-matrix(footprintCoverage,gidLn,1)
cmFootprint<-signalCounterMatrixTri(nSites=ZP,focalRange=fr,triTransform=gm$gid,signal=sig,sqDF=gm,mask=gm$exclude)
cmFootprintr<-signalCounterMatrixTri(nSites=ZPr,focalRange=fr,triTransform=gm$gid,signal=sig,sqDF=gm,mask=gm$exclude)

sig<-matrix(ojrWTrep$forClobber,gidLn,1)
cmOjrFor<-signalCounterMatrixTri(nSites=ZP,focalRange=fr,triTransform=gm$gid,signal=sig,sqDF=gm,mask=gm$exclude)

sig<-matrix(ojrWTrep$revClobber,gidLn,1)
cmOjrRev<-signalCounterMatrixTri(nSites=ZP,focalRange=fr,triTransform=gm$gid,signal=sig,sqDF=gm,mask=gm$exclude)


sig<-matrix(gm2$polCount,gidLn,1)
cmPolCount<-signalCounterMatrixTri(nSites=ZP,focalRange=fr,triTransform=gm$gid,signal=sig,sqDF=gm,mask=gm$exclude)

sig<-matrix(gm2$polCal,gidLn,1)
cmPolCall<-signalCounterMatrixTri(nSites=ZP,focalRange=fr,triTransform=gm$gid,signal=sig,sqDF=gm,mask=gm$exclude)

par(mar=c(4,5,4,5))
plot(fr,cmFootprint$Nsignal,type="n",axes=F,xlab="",ylab="")
axis(1)
polygon(c(fr,max(fr),min(fr)),c((cmFootprint$Nsignal),0,0),col="grey",border=NA)
#abline(v=0,col="lightblue")
par(new=T)
plot(fr,(cmOjrFor$NormRate+cmOjrRev$NormRate)/2,col="darkblue",type="n",xlab="",ylab="",axes=F)
points(fr,(cmOjrFor$NormRate+cmOjrRev$NormRate)/2,col="darkblue",lwd=lw,type="l")
axis(2,col="darkblue",col.axis="darkblue")
par(new=T)
ssPol<-smooth.spline(fr,cmPolCount$NormRate/cmPolCall$NormRate,df=35)
plot(ssPol,col="red",lwd=lw,type="l",axes=F,xlab="",ylab="")
axis(4,col="red",col.axis="red")
mtext("DNase footprint mid (nt)",side=1,line=2.5,col="black")
mtext("Okazaki junction rate (per nt)",side=2,line=3,col="darkblue")
mtext("Polymorphism rate",side=4,line=3,col="red")

dev.off()