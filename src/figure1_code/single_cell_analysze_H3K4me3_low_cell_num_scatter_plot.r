library(utils)
library(chromVAR)
library(motifmatchr)
library(Matrix)
library(SummarizedExperiment)
library(BiocParallel)
set.seed(2017)

register(MulticoreParam(16))
register(MulticoreParam(16, progressbar = TRUE))

##1-Figure 1b, 2-Extended Data Fig 4a, 3-Extended Data Fig 4b 4-Extended Data Fig 4c 5-Extended Data Fig 4d

kkk<-5

peaks <- getPeaks('./data/output/Figure1/KZ730-GA6909_0_0_mapq10_noDup-W200-G200-E.00000001.scoreisland',sort_peaks = TRUE)
peaks <-unique(peaks)
peaks<-sort(peaks)

if(kkk==1)
{
	aa2<-read.table("./data/input/Figure1/scChIC_3000a_3T3_readcount_at_3T3_peak.txt")
}  
if(kkk==2)
{
	aa2<-read.table("./data/input/Figure1/scChIC_3000b_3T3_readcount_at_3T3_peak.txt")

} 
if(kkk==3)
{
	aa2<-read.table("./data/input/Figure1/scChIC_1000b_readcount_at_3T3_peak.txt")

} 
if(kkk==4)
{
	aa2<-read.table("./data/input/Figure1/scChIC_300b_3T3_readcount_at_3T3_peak.txt")

} 
if(kkk==5)
{
	aa2<-read.table("./data/input/Figure1/scChIC_100b_3T3_readcount_at_3T3_peak.txt")

}                            


                             
aa1<-read.table("./data/input/Figure1/bulk_3T3_readcount_at_3T3_peak.txt")                                                              
    
                                                          
                                                          
bulk_cpm<-aa1[,1]*1000000/sum(aa1) 
sc_cpm<-aa2[,1]*1000000/sum(aa2)                             

xdata<-data.frame(log2(bulk_cpm[sc_cpm>0]+1),log2(sc_cpm[sc_cpm>0]+1))
acor<-cor(xdata, method="pearson")

#################################              
if(kkk==1)
{
	pdf("./Figures/Figure1/scatter_3T3_3000cell_at_H3K4me3_peaks_Fig1.pdf")
	smoothScatter(xdata, nbin = 300,
              colramp = colorRampPalette(c("white", "blue")),
              nrpoints = 300, pch = "", cex = 1.5, col = "black", xlab="Bulk Cell ChIP-seq (log2 CPM)", ylab ="3000 cells scChIC-seq (log2 CPM)",main = "", xlim=c(3.5,12), ylim=c(3.5,12),,cex.main=1.5, cex.lab=1.5,cex.axis=1.5)
	text(9,10,labels=paste0("Correlation = ", toString(format(round(acor[1,2], 2), nsmall = 2))), cex=1.5)	
	dev.off() 
}
if(kkk==2)
{
	pdf("./Figures/Figure1/scatter_3T3_3000cell_at_H3K4me3_peaks_Extfig.pdf")
	smoothScatter(xdata, nbin = 300,
              colramp = colorRampPalette(c("white", "blue")),
              nrpoints = 300, pch = "", cex = 1.5, col = "black", xlab="Bulk Cell ChIP-seq (log2 CPM)", ylab ="3000 cells scChIC-seq (log2 CPM)",main = "", xlim=c(3.5,12), ylim=c(3.5,12),,cex.main=1.5, cex.lab=1.5,cex.axis=1.5)
	text(9,10,labels=paste0("Correlation = ", toString(format(round(acor[1,2], 2), nsmall = 2))), cex=1.5)
	dev.off() 
}
if(kkk==3)
{
	pdf("./Figures/Figure1/scatter_3T3_1000cell_at_H3K4me3_peaks_Extfig.pdf")
	smoothScatter(xdata, nbin = 300,
              colramp = colorRampPalette(c("white", "blue")),
              nrpoints = 300, pch = "", cex = 1.5, col = "black", xlab="Bulk Cell ChIP-seq (log2 CPM)", ylab ="1000 cells scChIC-seq (log2 CPM)",main = "", xlim=c(3.5,12), ylim=c(3.5,12),,cex.main=1.5, cex.lab=1.5,cex.axis=1.5)
	text(9,10,labels=paste0("Correlation = ", toString(format(round(acor[1,2], 2), nsmall = 2))), cex=1.5)
	dev.off() 
}
if(kkk==4)
{
	pdf("./Figures/Figure1/scatter_3T3_300cell_at_H3K4me3_peaks_Extfig.pdf")
	smoothScatter(xdata, nbin = 300,
              colramp = colorRampPalette(c("white", "blue")),
              nrpoints = 300, pch = "", cex = 1.5, col = "black", xlab="Bulk Cell ChIP-seq (log2 CPM)", ylab ="300 cells scChIC-seq (log2 CPM)",main = "", xlim=c(3.5,12), ylim=c(3.5,12),,cex.main=1.5, cex.lab=1.5,cex.axis=1.5)
	text(9,10,labels=paste0("Correlation = ", toString(format(round(acor[1,2], 2), nsmall = 2))), cex=1.5)
	dev.off() 
}
if(kkk==5)
{
	#pdf("./Figures/Figure1/scatter_3T3_100cell_at_H3K4me3_peaks_Extfig.pdf")
	smoothScatter(xdata, nbin = 300,
              colramp = colorRampPalette(c("white", "blue")),
              nrpoints = 300, pch = "", cex = 1.5, col = "black", xlab="Bulk Cell ChIP-seq (log2 CPM)", ylab ="100 cells scChIC-seq (log2 CPM)",main = "", xlim=c(3.5,12), ylim=c(3.5,12),,cex.main=1.5, cex.lab=1.5,cex.axis=1.5)
	text(9,10,labels=paste0("Correlation = ", toString(format(round(acor[1,2], 2), nsmall = 2))), cex=1.5)
	dev.off() 
}
