library(utils)
library(chromVAR)
library(motifmatchr)
library(Matrix)
library(SummarizedExperiment)
library(BiocParallel)
set.seed(2017)

register(MulticoreParam(16))
register(MulticoreParam(16, progressbar = TRUE))

bulk0<-0

peaks <- getPeaks('./data/input/Figure1/KZ730-GA6909_0_0_mapq10_noDup-W200-G200-E.00000001.scoreisland',sort_peaks = TRUE)

peaks <-unique(peaks)
peaks<-sort(peaks)


if(bulk0==1)
{
	bulk_bedfile2<-paste0("./data/input/Figure1/KZ730-GA6909_0_0_mapq10_noDup.bed")
	bulk_counts<- getCounts(bulk_bedfile2, 
                             peaks, 
                             paired =  FALSE, 
                             by_rg = FALSE, 
                             format = "bed", 
                             colData = DataFrame(celltype ="bulk" ))
    xx0<-as.matrix(assays(bulk_counts)$counts)                     
    write.table(xx0,"./data/input/Figure1/bulk_3T3_readcount_at_3T3_peak.txt")                                                              
    
}            

bedfile3000a<-paste0("./GSE105012/","GSM3475877_KZ1252_GB0640_0_0_mapq10_noDup.bed")
bedfile3000b<-paste0("./GSE105012/","GSM3475847_GB0451_0_0_mapq10_noDup.bed")
bedfile1000b<-paste0("./GSE105012/","GSM3475849_GB0475_0_0_mapq10_noDup.bed")
bedfile300b<-paste0("./GSE105012/","GSM3475851_GB0476_0_0_mapq10_noDup.bed")
bedfile100b<-paste0("./GSE105012/","GSM3475853_GB0477_0_0_mapq10_noDup.bed")
                          


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
                                                                                                     
                                                              
                                                                                       
                                                 
         
xx1<-as.matrix(assays(fragment_counts_3000a)$counts)
xx2<-as.matrix(assays(fragment_counts_3000b)$counts)
xx3<-as.matrix(assays(fragment_counts_1000b)$counts)
xx4<-as.matrix(assays(fragment_counts_300b)$counts)
xx5<-as.matrix(assays(fragment_counts_100b)$counts)



write.table(xx1,"./data/temp/Figure1/scChIC_3000a_3T3_readcount_at_3T3_peak.txt")                                                              
write.table(xx2,"./data/temp/Figure1/scChIC_3000b_3T3_readcount_at_3T3_peak.txt")                                                              
write.table(xx3,"./data/temp/Figure1/scChIC_1000b_readcount_at_3T3_peak.txt")                                                              
write.table(xx4,"./data/temp/Figure1/scChIC_300b_3T3_readcount_at_3T3_peak.txt")                                                              
write.table(xx5,"./data/temp/Figure1/scChIC_100b_3T3_readcount_at_3T3_peak.txt")

