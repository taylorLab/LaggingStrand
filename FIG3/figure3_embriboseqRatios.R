###
# Load data from figuresDataLoader.R
strandRatios<-read.table(file="strandRatiosMasked.20140920.Rtable",header=T)
origins<-c(23920,67664,99593,113717,204507,228845,337352,374939,417363,442698,455372,540608,613043,684011,745349)

###

####
# Dependencies
library("Cairo")

####
# Pallet
delta<-rgb(10,6,0,maxColorValue=10) # delta
alpha<-rgb(8,0,0,maxColorValue=10)
epsilon<-rgb(0,6,10,maxColorValue=10)
wt<-rgb(5,5,5,maxColorValue=10)
wtwt<-rgb(2,2,2,maxColorValue=10)
epsilonSuper<-rgb(8,5,0,maxColorValue=10)
####

###
# chromosome 10 plots
ch10min<-min(which(gm$chr==10))
ch10max<-max(which(gm$chr==10))
ch10Windows<-which(strandRatios$gid>ch10min & strandRatios$gid<ch10max)
cXpos<-gm$chromPos[strandRatios$gid[ch10Windows]]
chr10length<-ch10max-ch10min
###

x11(width=10,height=7)
#Cairo(file="chrXOriginsAndStrandRatiosMasked.pdf",type="pdf",height=175,width=184,units="mm",bg="white",canvas="white",pointsize=10)
lw<-1
layout(matrix(c(5,1,1,2,2,3,4,6,6,7,8,8,9,10),14,1))
par(mar=c(0,6,0,6))
ctg<-which(gm$chr==10)
cxMax<-max(gm$chromPos[ctg])
#origins<-gm$chromPos[allOrigins[which(ars$V1=="chrX")]]
##
# Origins and external data
okazakiRatio<-log2((ojrWTrep$forDDbg[ctg]+1)/(ojrWTrep$revDDbg[ctg]+1))
sm<-smooth.spline(gm$chromPos[ctg],okazakiRatio,df=80)
predOJline<-predict(sm,cXpos)
plot(predOJline$x,predOJline$y,type="n",axes=F,xlab="",ylab="",xlim=c(1,cxMax))
abline(v=origins,col=rgb(10,0,0,5,maxColorValue=10),lwd=lw,lty=3)
abline(h=0,lwd=1,col="black")
points(predOJline$x,predOJline$y,type="l",lwd=lw,col="darkblue")
axis(2,at=c(-3,0,3),lwd=1)
mtext("Okazaki seq.",side=2,line=3)
axis(4,at=c(3,-3),labels=c("+ve","-ve"),lwd=1)
mtext("Lagging strand",side=4,line=3)

# delta versus epsilon
sm<-smooth.spline(cXpos,strandRatios$run25[ch10Windows],df=80)
plot(sm,type="n",axes=F,xlab="",ylab="",xlim=c(1,cxMax))
abline(h=0,lwd=1,col="black")
abline(v=origins,col=rgb(10,0,0,5,maxColorValue=10),lwd=lw,lty=3)
points(sm,type="l",lwd=lw,col=delta)
axis(2,at=c(-0.5,0,0.5),lwd=1,col=delta,col.axis=delta)
sm2<-smooth.spline(cXpos,strandRatios$run28[ch10Windows],df=80)
par(new=T)
plot(sm2,type="n",axes=F,xlab="",ylab="",xlim=c(1,cxMax))
points(sm2,type="l",lwd=lw,col=epsilon)
axis(4,at=c(-1.5,0,1.5),lwd=1,col=epsilon,col.axis=epsilon)
mtext(expression(paste(delta,"*",sep="")),side=2,line=3,col=delta,las=1,cex=1.8)
mtext(expression(paste(epsilon,"*",sep="")),side=4,line=3,col=epsilon,las=1,cex=1.8)

##
# Rag data
ctg<-which(gm$chr==10)
plot(gm$chromPos[ctg],gm$repTime[ctg],type="n",xlim=c(1,cxMax),ylim=c(20,55),axes=F,xlab="",ylab="")

points(gm$chromPos[ctg],gm$repTime[ctg],type="l",lwd=lw,col="black")
abline(v=origins,col=rgb(10,0,0,5,maxColorValue=10),lwd=lw,lty=3)
axis(2,at=c(20,50),labels=c("Late","Early"))
mtext("Replication time",side=2,line=3)
axis(1,at=c(1,chr10length),labels=F,lwd=1)
axis(1,at=c(1,2e5,4e5,6e5),labels=T,lwd=1)
mtext("Position on chrX (nt)",side=1,line=3)

plot(1,type="n",axes=F,main="",xlab="",ylab="")
plot(1,type="n",axes=F,main="",xlab="",ylab="")



##
# Alpha versus WT non-rescaled
sm<-smooth.spline(cXpos,strandRatios$run36[ch10Windows],df=80)
plot(sm,type="n",axes=F,xlab="",ylab="",xlim=c(1,cxMax))
abline(h=0,lwd=1,col="black")
#abline(v=origins,col="red",lwd=2,lty=3)
#axis(1,at=c(0,8)*1e5,labels=F,lwd=3)
#axis(1,at=c(0,2e5,4e5,6e5),labels=T,lwd=3)
#text(0.5,1.5,expression(epsilon),col=delta,cex=1.8)
sm2<-smooth.spline(cXpos,strandRatios$run35[ch10Windows],df=80)
points(sm2,type="l",lwd=lw,col=wt)
points(sm,type="l",lwd=lw,col=alpha)
#abline(v=origins,col=rgb(10,0,0,5,maxColorValue=10),lwd=lw,lty=3)
axis(1,at=c(1,chr10length),labels=F,lwd=1)
axis(1,at=c(1,2e5,4e5,6e5),labels=T,lwd=1)
axis(2,at=c(-0.4,0,0.4),lwd=1,col=alpha,col.axis=alpha)
axis(4,at=c(-0.4,0,0.4),lwd=1,col=wt,col.axis=wt)

mtext(expression(paste(alpha,"*",sep="")),side=2,line=3,col=alpha,las=1,cex=1.8)
mtext("POL",side=4,line=3,col=wt,las=1,cex=1)


mtext("Position on chrX (nt)",side=1,line=3)

plot(1,type="n",axes=F,main="",xlab="",ylab="")


sm<-smooth.spline(cXpos,strandRatios$run38[ch10Windows],df=80)
plot(sm,type="n",axes=F,xlab="",ylab="",xlim=c(1,cxMax))
abline(h=0,lwd=1,col="black")
sm2<-smooth.spline(cXpos,strandRatios$run37[ch10Windows],df=80)
points(sm2,type="l",lwd=lw,col=wt)
points(sm,type="l",lwd=lw,col=alpha)
axis(1,at=c(1,chr10length),labels=F,lwd=1)
axis(1,at=c(1,2e5,4e5,6e5),labels=T,lwd=1)
axis(2,at=c(-0.4,0,0.4),lwd=1,col=alpha,col.axis=alpha)
axis(4,at=c(-0.4,0,0.4),lwd=1,col=wt,col.axis=wt)

mtext(expression(paste(alpha,"*",sep="")),side=2,line=3,col=alpha,las=1,cex=1.8)
mtext("POL",side=4,line=3,col=wt,las=1,cex=1)


mtext("Position on chrX (nt)",side=1,line=3)


plot(1,type="n",axes=F,main="",xlab="",ylab="")

dev.off()

