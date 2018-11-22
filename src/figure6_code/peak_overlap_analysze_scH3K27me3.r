library(utils)
library(chromVAR)
library(motifmatchr)
library(Matrix)
library(SummarizedExperiment)
library(BiocParallel)
library(devtools)
library(Biobase)
library(preprocessCore)

set.seed(2017)

register(MulticoreParam(16))
register(MulticoreParam(16, progressbar = TRUE))


#peaks <- getPeaks('/data/kuw/biocore/wlku/pipeline/run1061_KZ1644_H3K27me3/sel_cell_0_2_SDS/SRR3503786_0_0_mapq10_noDup-W200-G600-E0.01.scoreisland',sort_peaks = TRUE)
#peaks <- getPeaks('/data/kuw/biocore/wlku/pipeline/human_WBC_H3K27me3/selbed/combined_peak.txt',sort_peaks = TRUE)
#peaks <- getPeaks('/data/kuw/biocore/wlku/pipeline/human_WBC_H3K27me3/selbed/combined_peak-E.0001.txt',sort_peaks = TRUE)
peaks <- getPeaks('./data/input/Figure6/combined_peaks_simplify_36224_2.txt' ,sort_peaks = TRUE)

peaks <-unique(peaks)
peaks<-sort(peaks)



#scpeaks <- getPeaks('./data/input/Figure6/combined_preci_012-W200-G600-E1000.scoreisland',sort_peaks = TRUE)
#scpeaks <- getPeaks('./data/temp/Figure6/combined_preci_012-W200-G600-E1000.scoreisland',sort_peaks = TRUE)
scpeaks <- getPeaks('./data/temp/Figure6/H3K27me3_combined_106-W200-G600-E1000.scoreisland',sort_peaks = TRUE)

#scpeaks <- getPeaks('/data/kuw/biocore/wlku/pipeline/run1070_KZ1655_KZ1656/temp_selbed/sel_st300000_bed/combined_st300000-W200-G600-E500.scoreisland',sort_peaks = TRUE)
#scpeaks <- getPeaks('/data/kuw/biocore/wlku/pipeline/run1070_KZ1655_KZ1656/temp_selbed/run1061_run1070_combined-W200-G600-E500.scoreisland',sort_peaks = TRUE)
scpeaks <-unique(scpeaks)
scpeaks<-sort(scpeaks)

overlap_pool<-findOverlaps(peaks, scpeaks)


xd_pool<-dim(as.matrix((queryHits(overlap_pool))))

pd1<-dim(as.matrix(width(scpeaks)))
pd2<-dim(as.matrix(width(peaks)))

x<- "Bulk wbc cell, number of H3K27me3 peaks = "
print(x)
print(pd2[1])
x<- "Pooled cellss number of H3K27me3 peaks = "
print(x)
print(pd1[1])
x<- "Overlap_peaks = "
print(x)
print(xd_pool[1])