library(zoo)
#library(ff)
#library(ffbase)
source("signalCounterMatrixTri.R")

nucsff<-read.table.ffdf(file="jiangPugh2009dyadV64_annotated.gidTab",header=T)
##
# Data loading
nucs<-read.table("jiangPugh2009dyadV64_annotated.gidTab",header=T)
reb<-read.table(file="gidReb1Unified.tab",header=T)
#ojrWTrep<-read.table("GSM835651_wt_replicate.V59toV64.cov.gid.ojr.Rtable",header=T)
ojrWTrep<-read.table(file="ojrWTrep.with2kAverages.table",header=T)
gm<-read.table("gm.20140919.Rtable",header=T)
dnase<-read.table("Hesselberth_2009_DNaseI_hypersensitive_sites_V64.bed.gidTab",header=F)
#triShuff<-read.table.ffdf("triPosShuffled.20140928.Rtable",header=T)
triShuff<-triCompShuffle
##
# Temporary
backgroundWindow<-2001
#ojrWTrep$forDD2kbg<-rollapply(ojrWTrep$forCovDD,width=backgroundWindow,by=1,FUN=mean,align="center",na.rm=FALSE,fill=NA)
#ojrWTrep$revDD2kbg<-rollapply(ojrWTrep$revCovDD,width=backgroundWindow,by=1,FUN=mean,align="center",na.rm=FALSE,fill=NA)
#write.table(ojrWTrep,file="ojrWTrep.with2kAverages.table")
##

lw<-1
lineSubr<-rgb(0,0,0,maxColorValue=10)
library(Cairo)
Cairo(file="fig1_OJRvsPol.pdf",type="pdf",height=160,width=180,units="mm",bg="white",canvas="white",pointsize=10)


#x11(width=10,height=6)

layout(matrix(c(1,1,3,1,1,4,2,2,5,2,2,6),3,4,byrow=F))
#layout.show(12)

# Nucleosome pol rate
par(mar=c(5,5,4,5))
source("FIG1_SYMSHAPE/nucPol.R")
mtext("a",side=3,adj=0,line=1,cex=1.5)

# Reb pol rate
par(mar=c(5,5,4,5))
source("FIG1_SYMSHAPE/rebPol.R")
mtext("b",side=3,adj=0,line=1,cex=1.5)


# Nuc triComp matching
par(mar=c(5,5,3,1))
source("FIG1_SYMSHAPE/nucTriComp.R")
mtext("c",side=3,adj=0,line=1,cex=1.5)


# Reb triComp matching
par(mar=c(5,5,3,1))
source("FIG1_SYMSHAPE/rebTriComp.R")
mtext("d",side=3,adj=0,line=1,cex=1.5)


# Nucleosome zoo rate
par(mar=c(5,2,3,5))
source("FIG1_SYMSHAPE/nucZoo.R")
mtext("e",side=3,adj=0,line=1,cex=1.5)


# Reb zoo rate
par(mar=c(5,2,3,5))
source("FIG1_SYMSHAPE/rebZoo.R")
mtext("f",side=3,adj=0,line=1,cex=1.5)


dev.off()
#dev.copy2pdf(file="nucAndRebHumps.pdf")

