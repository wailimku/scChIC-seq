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

#======================
##bulk_bedfile<-read.table("/data/kuw/biocore/wlku/pipeline/run1061_KZ1644_H3K27me3/sel_cell_0_2_SDS/combined_IB2_allsc.bed")
##bulk_bedfile2<-read.table("/data/kuw/biocore/wlku/pipeline/run1061_KZ1644_H3K27me3/sel_cell_wotrim/combined_allsc_wotrim.bed")

#bulk_bedfile<-read.table("/data/kuw/biocore/wlku/pipeline/human_WBC_H3K27me3/selbed/encode_sel18wbc_list.txt")
#bulk_bedfile2<-paste0("/data/kuw/biocore/wlku/pipeline/human_WBC_H3K27me3/selbed/", bulk_bedfile[,1])

#nx<-c("bulk")
#for (i in 1:18)
#{
#	nx[i]<-"bulk"
#}

#bulk_counts<- getCounts(bulk_bedfile2, 
#                             peaks, 
#                             paired =  FALSE, 
#                             by_rg = FALSE, 
#                             format = "bed", 
#                             colData = DataFrame(celltype =nx ))
#                             
                             
                             
#peaks_id <- filterPeaks(bulk_counts ,min_fragments_per_peak = 1, non_overlapping = TRUE,
#  ix_return = TRUE) 
#bulk_counts2 <- filterPeaks(bulk_counts ,min_fragments_per_peak = 1, non_overlapping = TRUE,
#  ix_return = FALSE)                               

#write.table(as.matrix(assays(bulk_counts2)$counts),'./data/input/Figure6/bulk_cells_18_chipseq_counts_at_WBCpeaks.txt');



aa2<-read.table("./data/input/Figure6/bulk_cells_18_chipseq_counts_at_WBCpeaks.txt");



sc_bedfile<-read.table("./data/temp/Figure6/sel_filterlist_106.txt")
sc_bedfile2<-paste0("./data/temp/Figure6/filtered_bed/", sc_bedfile[,2])
#sc_bedfile<-read.table("/data/kuw/biocore/wlku/pipeline/run1070_KZ1655_KZ1656/temp_selbed/sel_preci_012_bed/sel_106_bedlist.txt")
#sc_bedfile2<-paste0("/data/kuw/biocore/wlku/pipeline/run1070_KZ1655_KZ1656/temp_selbed/sel_preci_012_bed/", sc_bedfile[,2])


nx<-c("sc")
for (i in 1:106)
{
	nx[i]<-"sc"
}
sincell_counts<- getCounts(sc_bedfile2, 
                             peaks, 
                             paired =  FALSE, 
                             by_rg = FALSE, 
                             format = "bed", 
                             colData = DataFrame(celltype =nx ))
                             
                             
sc_counts<-as.matrix(assays(sincell_counts)$counts)
aa<-sc_counts
aa[sc_counts>1]<-1                           
sc_counts<-aa                             
                             
sc_counts<-as.matrix(assays(sincell_counts)$counts)
sc_cpm<-log2(rowSums(sc_counts)*1000000/(sum(sc_counts))+1)
sc_cpm<-rowMeans(log2(sc_counts*1000000/(colSums(sc_counts)+1)+1))


yy <-colSums(sc_counts)/colData(sincell_counts)$depth
xx3<-sc_counts[,yy>0&colData(sincell_counts)$depth>1000] 

sc_cpm<-rowMeans(log2(xx3*1000000/(colSums(xx3)+0.0001)+1))


acounts<- as.matrix(aa2)
bulk_cpm<-rowMeans(acounts[,c(1:18)]*1000000/(colSums(acounts[,c(1:18)])))

xdata<-data.frame(sc_cpm,log2(bulk_cpm+1))
norm_data = normalize.quantiles(as.matrix(xdata[bulk_cpm>3&sc_cpm>0,]))
acor<-cor(norm_data, method="pearson")



xdata<-data.frame(sc_cpm,bulk_cpm)


lib_read<-colSums(as.matrix(assays(sincell_counts)$counts))
lib_total<-sc_bedfile[,1]
preci<- lib_read/lib_total          
aa<-as.matrix(assays(sincell_counts)$counts)
sc_sen<-as.matrix(assays(sincell_counts)$counts)
sc_sen[aa>0]<-1
sensi<-colSums(sc_sen)*100/36224


                             
                             
cc_counts<-as.matrix(assays(sincell_counts)$counts)                             
lib_read<-colSums(as.matrix(assays(sincell_counts)$counts))
lib_total<-sc_bedfile[,1]
preci<- lib_read/lib_total  

write.table(cc_counts,'./data/temp/Figure6/bulk_wbc_rc_36224_106_mat_2.txt')
write.table(as.matrix(colnames(cc_counts)),'./data/temp/Figure6/sc_wbc_file_name_2.txt')
#write.table(data.frame(seqnames(peaks),ranges(peaks)),'./data/temp/Figure6/combind_peaks_simplify_36224_2.txt')          
                            
                            
pdf("./Figures/Figure6/scatter_bulk_pooledsc_at_H3K27me3peaks_3.pdf")
       
smoothScatter(norm_data, nbin = 300,
              colramp = colorRampPalette(c("white", "blue")),
              nrpoints = 300, pch = "", cex = 1, col = "black", xlab="pool scWBC log2(CPM)", ylab ="Bulk WBC log2(CPM)",main = "", xlim=c(1.5,15), ylim=c(1.5,15),cex.main=1.5, cex.lab=1.5,cex.axis=1.5)
                  
text(9,10,labels=paste0("Correlation = ", toString(format(round(acor[1,2], 2), nsmall = 2))), cex=1.5)
dev.off()                             
                             
#======================
