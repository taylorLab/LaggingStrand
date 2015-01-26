# Martin Taylor
# (c) 2014
# Efficient signal counter for defined range around list of focal sites.
# Modified version allowing for a sites transfrom matrix to be provided.
# With this we can use perfectly compositionally matched (to abitrary k-mer)
# matricies of sites.

signalCounterMatrixTri<-function(nSites,focalRange,triTransform,signal,sqDF,mask){
 # nSites = GIDs of focal sites.
 # focalRange = the range around focal sites to calculate
 # signal = the counts/signal vector
 # sqDF = reference to the sq dataframe
 # mask = masking vector (1=masked)
 frq<-data.frame(a=c(rep(0,length(focalRange))),t=c(rep(0,length(focalRange))),c=c(rep(0,length(focalRange))),g=c(rep(0,length(focalRange))))
 frq$an<-0
 frq$tn<-0
 frq$cn<-0
 frq$gn<-0
 frq$focal<-NA
 maxGid<-max(sqDF$gid)
 forOrder<-c("A","C","G","T")
 theOrder<-forOrder
 cols<-dim(signal)[2]
 for(i in 1:length(focalRange)){
  ii<-focalRange[i]
  frq$focal[i]=ii
  qsite<-nSites+ii
  ###
  qsite[which(qsite<1)]=1
  qsite[which(qsite>maxGid)]=maxGid
  ####
  # The triCompShuffle is brought in here (can just provide gid for no transfrom)
  qsite<-triTransform[qsite]
  ####
  # honour exclusions
  qsite=qsite[which(mask[qsite]==0)]
  ##
  ntA<-qsite[which(sqDF$forNuc[qsite]==theOrder[1])]
  frq$an[i]=sum(!is.na(signal[ntA,c(1:cols)]))
  frq$a[i]=sum(signal[ntA,c(1:cols)],na.rm=T)
  
  ntT<-qsite[which(sqDF$forNuc[qsite]==theOrder[4])]
  frq$tn[i]=sum(!is.na(signal[ntT,c(1:cols)]))
  frq$t[i]=sum(signal[ntT,c(1:cols)],na.rm=T)
  
  ntC<-qsite[which(sqDF$forNuc[qsite]==theOrder[2])]
  frq$cn[i]=sum(!is.na(signal[ntC,c(1:cols)]))
  frq$c[i]=sum(signal[ntC,c(1:cols)],na.rm=T)
  
  ntG<-qsite[which(sqDF$forNuc[qsite]==theOrder[3])]
  frq$gn[i]=sum(!is.na(signal[ntG,c(1:cols)]))
  frq$g[i]=sum(signal[ntG,c(1:cols)],na.rm=T)
 }
 frq$Nsignal<-frq$a+frq$t+frq$c+frq$g
 frq$Nrate<-frq$Nsignal/(frq$an+frq$tn+frq$cn+frq$gn)
 frq$NormRate<-((frq$a/frq$an)+(frq$t/frq$tn)+(frq$c/frq$cn)+(frq$g/frq$gn))/rep(4,length(frq$a))
 return(frq)
}
