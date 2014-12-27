#####
# Function

library(zoo)
peakFinder<-function(series,span=20,target=2){
 modSpan=((span%/%2)*2)+1
 midSpan=(span%/%2)+1
 rxz <- rollapply(series, modSpan, align="center", function(x) which.max(x)==midSpan)
 irxz<-index(rxz)[coredata(rxz)]
 rt<-irxz[which(rank(series[irxz],na.last=T)>length(irxz)-target)]
 rt+midSpan
}


#####



simlen<-1e9
gelIntensityThresh<-0
#promega1kbLadder<-c(10000,8000,6000,5000,4000,3000,2000,1500,1000)
promega1kbLadder<-c(10000,8000,6000,5000,4000,3000,2000,1500)



ladderYscale<-((max(c(mef$alpha)))/max(mef$m))/2
lanes<-c(2:length(mef))
samp<-c((dim(mef)[1]-100):dim(mef)[1])
sampVec<-c(mef[samp,3],mef[samp,4])
gelIntensityThresh<-as.numeric(quantile(sampVec,probs=c(0.9),na.rm=T))



mef$wt=mef$wt-gelIntensityThresh
mef$wt[which(mef$wt<0)]=0
mef$rnh=mef$rnh-gelIntensityThresh
mef$rnh[which(mef$rnh<0)]=0
mef$alpha=mef$alpha-gelIntensityThresh
mef$alpha[which(mef$alpha<0)]=0

x11()
layout(matrix(c(1,2,3,4),2,2))

#####
# Color pallet
colMarker<-"grey"
colPP<-"#1924b1"
colPM<-"#65a5d1"
colMM<-"#ff9000"
colMM2<-"#ffbf00"
colCUT<-"#00a876"
#####



#####################################################
#######
# ladder calibration
thrLadderVec<-mef$m
pvpSize<-promega1kbLadder
# pvp is the gel position of the peaks
pvpIndex<-peakFinder(thrLadderVec,span=30,target=length(pvpSize))
pvp<-mef$Position[pvpIndex]
# Fit linar model
x<-c(pvpSize[1:length(pvp)])
y<-c(pvp)
fit<-lm(y~log(x))
calibIntercept<-as.numeric(fit$coefficients[1])
calibSlope<-as.numeric(fit$coefficients[2])
mef$ntPos<-exp((mef$Position-calibIntercept)/calibSlope)
mef[which(mef$ntPos>30000),c(3,4,5)]=0
mef[which(mef$ntPos<1200),c(3,4,5)]=0

zeroMigrationPoint<-exp((0-calibIntercept)/calibSlope)
####
# Check calibration looks good
calcol<-"grey"
plot(pvpSize,pvp,xlab="Fragment length (nt)",ylab="Electrophoretic distance",ylim=c(0,74),xlim=c(0,zeroMigrationPoint),pch=20,cex=1)
 nw<-data.frame(x=sequence(zeroMigrationPoint))
 psl<-predict(fit,nw,interval="prediction")
 lines(sequence(zeroMigrationPoint),psl[,2],col="grey",lwd=3)
 lines(sequence(zeroMigrationPoint),psl[,3],col="grey",lwd=3)
 lines(sequence(zeroMigrationPoint),psl[,1],lwd=3)
 points(pvpSize[1:14],pvp[1:14],cex=2)

plot(mef$Position,mef$rnh,type="n",ylim=c(0,max(mef$wt)),xlab="Migration distance",ylab="Intensity")
 points(mef$Position,mef$m*ladderYscale,type="l",lwd=1,col="grey")
 points(mef$Position[pvpIndex],mef$m[pvpIndex]*ladderYscale,lwd=1,col="orange")
 points(mef$Position,mef$wt,type="l",lwd=3,col="black")
 points(mef$Position,mef$rnh,type="l",lwd=3,col="blue")
 points(mef$Position,mef$alpha,type="l",lwd=3,col="darkred")
#######
#####################################################


## extrapolate calibration curve for datapoints that have 
## not been observed (needed for simulation space)
emb<-embed(mef$ntPos,2)
embm<-emb[,2]-emb[,1]
x<-emb[,1]
y<-embm
fitSliceSize<-lm(y~x)
extendX<-sort(mef$ntPos)
newXlast<-max(mef$ntPos)
while(newXlast<zeroMigrationPoint){
 nxf<-data.frame(x=c(newXlast))
 nx<-predict(fitSliceSize,nxf,interval="prediction")
 newXlast<-nx[,1]+newXlast
 extendX<-c(extendX,newXlast)
}
newXfirst<-min(mef$ntPos)
while(newXfirst>1){
 nxf<-data.frame(x=c(newXfirst))
 nx<-predict(fitSliceSize,nxf,interval="prediction")
 newXfirst<-newXfirst-nx[,1]
 extendX<-c(newXfirst,extendX)
 #print(newXfirst)
}
eXt<-data.frame(x=c(extendX))
neX<-predict(fit,eXt,interval="prediction")
eX<-data.frame(nt=extendX,migration=neX[,1],templt=0)
##


sdf<-10
## WT
deint<-mef$wt/mef$ntPos
deint[is.na(deint)]<-0
Spline<-(smooth.spline(mef$ntPos,deint,df=sdf))
Hist<-predict(Spline,eX$nt,interval="prediction")
Hist$y[which(Hist$y<0)]<-0
Sum<-sum(Hist$y*eX$nt)
wtHist1G<-round((Hist$y/Sum)*simlen)
wtVals1G<-rep(eX$nt,wtHist1G)

## rnh201
deint<-mef$rnh/mef$ntPos
deint[is.na(deint)]<-0
Spline<-(smooth.spline(mef$ntPos,deint,df=sdf))
Hist<-predict(Spline,eX$nt,interval="prediction")
Hist$y[which(Hist$y<0)]<-0
Sum<-sum(Hist$y*eX$nt)
rnhHist1G<-round((Hist$y/Sum)*simlen)
rnhVals1G<-rep(eX$nt,rnhHist1G)

## Alpha
deint<-mef$alpha/mef$ntPos
deint[is.na(deint)]<-0
Spline<-(smooth.spline(mef$ntPos,deint,df=sdf))
Hist<-predict(Spline,eX$nt,interval="prediction")
Hist$y[which(Hist$y<0)]<-0
Sum<-sum(Hist$y*eX$nt)
alphaHist1G<-round((Hist$y/Sum)*simlen)
alphaVals1G<-rep(eX$nt,alphaHist1G)


### projecting in frag-space
ymax<-max(c(alphaHist1G))
plot(eX$nt,wtHist1G,ylab="Fragment count",xlab="Fragment length",ylim=c(0,ymax),type="n")
points(eX$nt,wtHist1G,type="l",col="black",lwd=3)
points(eX$nt,rnhHist1G,type="l",col="blue",lwd=3)
points(eX$nt,alphaHist1G,type="l",col="darkred",lwd=3)

genomeScaleFactor<-12e6/1e9
rnhRibo<-(sum(rnhHist1G)-sum(wtHist1G))*genomeScaleFactor
alphaRibo<-(sum(alphaHist1G)-sum(wtHist1G))*genomeScaleFactor


plot(c(1:length(eX$nt)),wtHist1G,ylab="Fragment count",xlab="Fragment length",ylim=c(0,ymax),type="n")
points(c(1:length(eX$nt)),wtHist1G,type="l",col="black",lwd=3)
points(c(1:length(eX$nt)),rnhHist1G,type="l",col="blue",lwd=3)
points(c(1:length(eX$nt)),alphaHist1G,type="l",col="darkred",lwd=3)

