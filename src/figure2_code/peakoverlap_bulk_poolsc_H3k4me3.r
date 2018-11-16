library(chromVAR)
library(motifmatchr)
library(Matrix)
library(SummarizedExperiment)
library(BiocParallel)
library(utils)

set.seed(2017)
register(MulticoreParam(8))
register(MulticoreParam(8, progressbar = TRUE))


peaks <- getPeaks("./data/output/Figure2/peakname3.txt", sort_peaks = TRUE)
scpeaks <- getPeaks("./data/output/Figure2/singlecell_read1-W200-G200-E.01.scoreisland", sort_peaks = TRUE)


overlap_pool<-findOverlaps(peaks, scpeaks)


xd_pool<-dim(as.matrix((queryHits(overlap_pool))))

pd1<-dim(as.matrix(width(scpeaks)))
pd2<-dim(as.matrix(width(peaks)))

	x<- "Bulk wbc cell, number of H3K4me3 peaks = "
	print(x)
	print(pd2[1])
	x<- "Pooled single cells scChIC-seq, number of H3K4me3 peaks = "
	print(x)
	print(pd1[1])
	x<- "Overlap_peaks = "
	print(x)
	print(xd_pool[1])