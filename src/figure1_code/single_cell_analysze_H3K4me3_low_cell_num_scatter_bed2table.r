library(utils)
library(chromVAR)
library(motifmatchr)
library(Matrix)
library(SummarizedExperiment)
library(BiocParallel)
set.seed(2017)

register(MulticoreParam(16))
register(MulticoreParam(16, progressbar = TRUE))

peaks <- getPeaks('./data/input/Figure1/KZ730-GA6909_0_0_mapq10_noDup-W200-G200-E.00000001.scoreisland',sort_peaks = TRUE)

peaks <-unique(peaks)
peaks<-sort(peaks)

bulk_bedfile2<-paste0("./data/input/Figure1/KZ730-GA6909_0_0_mapq10_noDup.bed")


bulk_counts<- getCounts(bulk_bedfile2, 
                             peaks, 
                             paired =  FALSE, 
                             by_rg = FALSE, 
                             format = "bed", 
                             colData = DataFrame(celltype ="bulk" ))
            

bedfile3000a<-paste0("./data/input/Figure1/","KZ1252_GB0640_0_0_mapq10_noDup.bed")
bedfile3000b<-paste0("./data/input/Figure1/","KZ1225_GB0451_0_0_mapq10_noDup.bed")
bedfile1000b<-paste0("./data/input/Figure1/","KZ1245-2_0_0_mapq10_noDup.bed")
bedfile300b<-paste0("./data/input/Figure1/","KZ1245-3_0_0_mapq10_noDup.bed")
bedfile100b<-paste0("./data/input/Figure1/","KZ1245-4_0_0_mapq10_noDup.bed")
                          


fragment_counts_3000a <- getCounts(bedfile3000a, 
                             peaks, 
                             paired =  FALSE, 
                             by_rg = FALSE, 
                             format = "bed", 
                             colData = DataFrame(celltype ="3T3"))
fragment_counts_3000b <- getCounts(bedfile3000b, 
                             peaks, 
                             paired =  FALSE, 
                             by_rg = FALSE, 
                             format = "bed", 
                             colData = DataFrame(celltype ="3T3"))
fragment_counts_1000b <- getCounts(bedfile1000b, 
                             peaks, 
                             paired =  FALSE, 
                             by_rg = FALSE, 
                             format = "bed", 
                             colData = DataFrame(celltype ="3T3"))                                                          
fragment_counts_300b <- getCounts(bedfile300b, 
                             peaks, 
                             paired =  FALSE, 
                             by_rg = FALSE, 
                             format = "bed", 
                             colData = DataFrame(celltype ="3T3"))
fragment_counts_100b <- getCounts(bedfile100b, 
                             peaks, 
                             paired =  FALSE, 
                             by_rg = FALSE, 
                             format = "bed", 
                             colData = DataFrame(celltype ="3T3"))
                                                                                                     
                                                              
                                                                                       
                                                 
xx0<-as.matrix(assays(bulk_counts)$counts)         
xx1<-as.matrix(assays(fragment_counts_3000a)$counts)
xx2<-as.matrix(assays(fragment_counts_3000b)$counts)
xx3<-as.matrix(assays(fragment_counts_1000b)$counts)
xx4<-as.matrix(assays(fragment_counts_300b)$counts)
xx5<-as.matrix(assays(fragment_counts_100b)$counts)



write.table(xx0,"./data/input/Figure1/bulk_3T3_readcount_at_3T3_peak.txt")                                                              
write.table(xx1,"./data/input/Figure1/scChIC_3000a_3T3_readcount_at_3T3_peak.txt")                                                              
write.table(xx2,"./data/input/Figure1/scChIC_3000b_3T3_readcount_at_3T3_peak.txt")                                                              
write.table(xx3,"./data/input/Figure1/scChIC_1000b_readcount_at_3T3_peak.txt")                                                              
write.table(xx4,"./data/input/Figure1/scChIC_300b_3T3_readcount_at_3T3_peak.txt")                                                              
write.table(xx5,"./data/input/Figure1/scChIC_100b_3T3_readcount_at_3T3_peak.txt")