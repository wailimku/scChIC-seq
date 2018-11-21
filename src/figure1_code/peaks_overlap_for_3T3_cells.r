#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(utils)
library(chromVAR)
library(motifmatchr)
library(Matrix)
library(SummarizedExperiment)
library(BiocParallel)
set.seed(2017)

register(MulticoreParam(16))
register(MulticoreParam(16, progressbar = TRUE))

############################
# kkk=1 for Fig 1d, kkk=2 for Extended Data Fig 1c, kkk=3 for Extended Data Fig 1d
############################
kkk<-args[1] #Plese change
############################
peaks <- getPeaks('./data/output/Figure1/KZ730-GA6909_0_0_mapq10_noDup-W200-G200-E.00000001.scoreisland',sort_peaks = TRUE)

peaks <-unique(peaks)
peaks<-sort(peaks)

chicpeaks <- getPeaks('./data/temp/Figure1/GSM3475877_KZ1252_GB0640_0_0_mapq10_noDup-W200-G200-E500.scoreisland',sort_peaks = TRUE)



chicpeaks3000 <- getPeaks('./data/temp/Figure1/GSM3475847_GB0451_0_0_mapq10_noDup-W200-G200-E100.scoreisland',sort_peaks = TRUE)
chicpeaks1000 <- getPeaks('./data/temp/Figure1/GSM3475849_GB0475_0_0_mapq10_noDup-W200-G200.scoreisland',sort_peaks = TRUE)
chicpeaks300 <- getPeaks('./data/temp/Figure1/GSM3475851_GB0476_0_0_mapq10_noDup-W200-G200.scoreisland',sort_peaks = TRUE)
chicpeaks100 <- getPeaks('./data/temp/Figure1/GSM3475853_GB0477_0_0_mapq10_noDup-W200-G200.scoreisland',sort_peaks = TRUE)



overlap3000<-findOverlaps(peaks, chicpeaks3000)
overlap1000<-findOverlaps(peaks, chicpeaks1000)
overlap300<-findOverlaps(peaks, chicpeaks300)
overlap100<-findOverlaps(peaks, chicpeaks100)

xd3000<-dim(as.matrix((queryHits(overlap3000))))
xd1000<-dim(as.matrix((queryHits(overlap1000))))
xd300<-dim(as.matrix((queryHits(overlap300))))
xd100<-dim(as.matrix((queryHits(overlap100))))

pd3000<-dim(as.matrix(width(chicpeaks3000)))
pd1000<-dim(as.matrix(width(chicpeaks1000)))
pd300<-dim(as.matrix(width(chicpeaks300)))
pd100<-dim(as.matrix(width(chicpeaks100)))


overlapchic<-findOverlaps(peaks, chicpeaks)
xdchic<-dim(as.matrix((queryHits(overlapchic))))

pd1<-dim(as.matrix(width(chicpeaks)))
pd2<-dim(as.matrix(width(peaks)))
if(kkk==1)
{
	x<- "Bulk cell, number of peaks = "
	print(x)
	print(pd2[1])
	x<- "3000 cells scChIC-seq, number of peaks = "
	print(x)
	print(pd1[1])
	x<- "Overlap_peaks = "
	print(x)
	print(xdchic[1])	
}


dda<-c(xd3000[1]/pd3000[1],xd1000[1]/pd1000[1],xd300[1]/pd300[1],xd100[1]/pd100[1])

if(kkk==2)
{
	pdf("./Figures/Figure1/barplot_precision_3T3_low_cell_num.pdf")
	barplot(dda, main="scChIC-seq for 3T3 ",ylab="Precison",xlab="Number of cells", names.arg=c("3000", "1000 ", "300","100 "), cex.names=1.5,,cex.main=1.5, cex.lab=1.5,cex.axis=1.5)
	dev.off() 
}

pda<-c(xd3000[1],xd1000[1],xd300[1],xd100[1])


if(kkk==3)
{
	pdf("./Figures/Figure1/barplot_senstivity_3T3_low_cell_num.pdf")
	barplot(pda, main="scChIC-seq for 3T3 ",ylab="Number of Peaks",xlab="Number of cells", names.arg=c("3000", "1000 ", "300","100 "), cex.names=1.5,,cex.main=1.5, cex.lab=1.5,cex.axis=1.5)
	dev.off() 
}




