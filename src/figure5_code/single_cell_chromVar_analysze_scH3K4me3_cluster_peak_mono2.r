library(chromVAR)
library(motifmatchr)
library(Matrix)
library(SummarizedExperiment)
library(BiocParallel)
library(utils)

set.seed(2017)
register(MulticoreParam(16))
register(MulticoreParam(16, progressbar = TRUE))


B1<- read.table('./data/input/Figure4/scimpute/scimpute_count.txt',header=TRUE,row.name="Gene")
#B1<- read.table('/data/kuw/biocore/wlku/Keji/scwbc_rc_52798_242_mat2.txt',header=TRUE,row.name="Gene")
peaks2 <- read.table("./data/input/Figure2/peakname3.txt", header=FALSE)
peaks <- getPeaks("./data/input/Figure2/peakname3.txt", sort_peaks = TRUE)
peaks<- resize(peaks, width = 3000, fix = "center")

#B1<-B1[,-112]
B1[,112]=NULL
sc3_cellclus <- read.table('./data/input/Figure4/cell_clus_sc3_12779_wimpute_nobatch_7a.txt')

B2<-B1[,sc3_cellclus==7]  #####need to change

#kmed_peak<- read.table("/data/kuw/biocore/wlku/Keji/wbc_285_read1_mapped/conet_het_result/conet_het_scrna_sch3k4_kmedoids_6_c6.txt",header=TRUE)
#kmed_peak<- read.table("/data/kuw/biocore/wlku/Keji/wbc_285_read1_mapped/conet_het_result/conet_het_scrna_sch3k4_rhigh_het_clow_het_mono3.txt",header=TRUE)
#kmed_peak<- read.table("/data/kuw/biocore/wlku/Keji/wbc_285_read1_mapped/conet_het_result/conet_het_scrna_sch3k4_rhigh_het_chigh_het2_mono.txt",header=TRUE)
kmed_peak<- read.table("./data/output/Figure5/conet_het_result/conet_het_scrna_sch3k4_rhigh_het_chigh_het2_mono.txt",header=TRUE)


uni_chr<-unique(peaks2[,1])
library(BSgenome.Hsapiens.UCSC.hg18)
genome <- BSgenome.Hsapiens.UCSC.hg18
m<-1
l<-{}
for (i in uni_chr)
{
	x<-grep(TRUE,peaks2[,1]==i)
	x1<-grep(TRUE,peaks2[x,2]<0 | peaks2[x,3]>seqlengths(genome)[m])
	m<-m+1
	if(is.nan(x1[1])==FALSE)
	{
		l<-append(l,x1+x[1]-1)
	}
}

peaks2<-peaks2[-l,]
peaks<-peaks[-l]
B2<-B2[-l,]


colData<-DataFrame(colnames(B2),colSums(as.matrix(B2)))
colnames(colData)<-c("cell","depth")
example_counts<-SummarizedExperiment(assays=list(counts=as.matrix(B2)),rowRanges=peaks,colData=colData)

library(BSgenome.Hsapiens.UCSC.hg18)
example_counts <- addGCBias(example_counts, genome = BSgenome.Hsapiens.UCSC.hg18)

counts_filtered <- filterSamples(example_counts, min_depth = 1000, 
                                 min_in_peaks = 0.15, shiny = FALSE)
                                 
                                 
ix <- filterSamples(example_counts, ix_return = TRUE, shiny = FALSE)
counts_filtered <- filterPeaks(counts_filtered, non_overlapping = TRUE)
                                 
id<-match(rownames(assays(counts_filtered)$counts),kmed_peak[,1])    
zz<-grep(TRUE,id>0)

counts_filtered2<-counts_filtered[zz,]
                             
library(chromVARmotifs)
data("human_pwms_v1") # human collection         
motifs<-human_pwms_v1                                 

motif_ix <- matchMotifs(motifs, counts_filtered2, genome = BSgenome.Hsapiens.UCSC.hg18)

bg <- getBackgroundPeaks(object = counts_filtered2)
dev <- computeDeviations(object = counts_filtered2, annotations = motif_ix,
                         background_peaks = bg)
variability <- computeVariability(dev)

#tsne_plots <- plotDeviationsTsne(dev, tsne_results,  shiny = FALSE)

aa<-rownames(dev)
#write.table(as.matrix(aa),file="/data/kuw/biocore/wlku/Keji/wbc_285_read1_mapped/conet_het_result/conet_het_scrna_sch3k4_rhigh_het_chigh_het2_motif_mono.txt",sep="\t",quote=FALSE,row.names=TRUE, na = "1000")
write.table(as.matrix(aa),file="./data/output/Figure5/conet_het_result/conet_het_scrna_sch3k4_rhigh_het_chigh_het2_motif_mono.txt",sep="\t",quote=FALSE,row.names=TRUE, na = "1000")

aa<-assays(dev)$z
#write.table(data.frame(aa),file="/data/kuw/biocore/wlku/Keji/wbc_285_read1_mapped/chromvar/conet_het_scrna_sch3k4_rhigh_het_chigh_het2_zscore_mono.txt",sep="\t",quote=FALSE,row.names=TRUE,na = "1000")
write.table(data.frame(aa),file="./data/output/Figure5/chromvar/conet_het_scrna_sch3k4_rhigh_het_chigh_het2_zscore_mono.txt",sep="\t",quote=FALSE,row.names=TRUE,na = "1000")

aa<-assays(dev)$deviations
#write.table(data.frame(aa),file="/data/kuw/biocore/wlku/Keji/wbc_285_read1_mapped/chromvar/conet_het_scrna_sch3k4_rhigh_het_chigh_het2_dev_mono.txt",sep="\t",quote=FALSE,row.names=TRUE,na = "1000")
#write.table(data.frame(variability ),file="/data/kuw/biocore/wlku/Keji/wbc_285_read1_mapped/chromvar/conet_het_scrna_sch3k4_rhigh_het_chigh_het2_var_mono.txt",sep="\t",quote=FALSE,row.names=TRUE,na = "1000")
write.table(data.frame(aa),file="./data/output/Figure5/chromvar/conet_het_scrna_sch3k4_rhigh_het_chigh_het2_dev_mono.txt",sep="\t",quote=FALSE,row.names=TRUE,na = "1000")
write.table(data.frame(variability ),file="./data/output/Figure5/chromvar/conet_het_scrna_sch3k4_rhigh_het_chigh_het2_var_mono.txt",sep="\t",quote=FALSE,row.names=TRUE,na = "1000")



#################################################################
#################################################################

library(chromVAR)
library(motifmatchr)
library(Matrix)
library(SummarizedExperiment)
library(BiocParallel)
library(utils)

set.seed(2017)
register(MulticoreParam(16))
register(MulticoreParam(16, progressbar = TRUE))


B1<- read.table('./data/input/Figure4/scimpute/scimpute_count.txt',header=TRUE,row.name="Gene")
#B1<- read.table('/data/kuw/biocore/wlku/Keji/scwbc_rc_52798_242_mat2.txt',header=TRUE,row.name="Gene")
peaks2 <- read.table("./data/input/Figure2/peakname3.txt", header=FALSE)
peaks <- getPeaks("./data/input/Figure2/peakname3.txt", sort_peaks = TRUE)
peaks<- resize(peaks, width = 3000, fix = "center")

#B1<-B1[,-112]
B1[,112]=NULL
sc3_cellclus <- read.table('./data/input/Figure4/cell_clus_sc3_12779_wimpute_nobatch_7a.txt')

B2<-B1[,sc3_cellclus==7]  #####need to change

#kmed_peak<- read.table("/data/kuw/biocore/wlku/Keji/wbc_285_read1_mapped/conet_het_result/conet_het_scrna_sch3k4_kmedoids_6_c6.txt",header=TRUE)
#kmed_peak<- read.table("/data/kuw/biocore/wlku/Keji/wbc_285_read1_mapped/conet_het_result/conet_het_scrna_sch3k4_rhigh_het_clow_het_mono3.txt",header=TRUE)
#kmed_peak<- read.table("/data/kuw/biocore/wlku/Keji/wbc_285_read1_mapped/conet_het_result/conet_het_scrna_sch3k4_rhigh_het_chigh_het2_mono.txt",header=TRUE)
kmed_peak<- read.table("./data/output/Figure5/conet_het_result/conet_het_scrna_sch3k4_rhigh_het_clow_het2_mono.txt",header=TRUE)


uni_chr<-unique(peaks2[,1])
library(BSgenome.Hsapiens.UCSC.hg18)
genome <- BSgenome.Hsapiens.UCSC.hg18
m<-1
l<-{}
for (i in uni_chr)
{
	x<-grep(TRUE,peaks2[,1]==i)
	x1<-grep(TRUE,peaks2[x,2]<0 | peaks2[x,3]>seqlengths(genome)[m])
	m<-m+1
	if(is.nan(x1[1])==FALSE)
	{
		l<-append(l,x1+x[1]-1)
	}
}

peaks2<-peaks2[-l,]
peaks<-peaks[-l]
B2<-B2[-l,]


colData<-DataFrame(colnames(B2),colSums(as.matrix(B2)))
colnames(colData)<-c("cell","depth")
example_counts<-SummarizedExperiment(assays=list(counts=as.matrix(B2)),rowRanges=peaks,colData=colData)

library(BSgenome.Hsapiens.UCSC.hg18)
example_counts <- addGCBias(example_counts, genome = BSgenome.Hsapiens.UCSC.hg18)

counts_filtered <- filterSamples(example_counts, min_depth = 1000, 
                                 min_in_peaks = 0.15, shiny = FALSE)
                                 
                                 
ix <- filterSamples(example_counts, ix_return = TRUE, shiny = FALSE)
counts_filtered <- filterPeaks(counts_filtered, non_overlapping = TRUE)
                                 
id<-match(rownames(assays(counts_filtered)$counts),kmed_peak[,1])    
zz<-grep(TRUE,id>0)

counts_filtered2<-counts_filtered[zz,]
                             
library(chromVARmotifs)
data("human_pwms_v1") # human collection         
motifs<-human_pwms_v1                                 

motif_ix <- matchMotifs(motifs, counts_filtered2, genome = BSgenome.Hsapiens.UCSC.hg18)

bg <- getBackgroundPeaks(object = counts_filtered2)
dev <- computeDeviations(object = counts_filtered2, annotations = motif_ix,
                         background_peaks = bg)
variability <- computeVariability(dev)

#tsne_plots <- plotDeviationsTsne(dev, tsne_results,  shiny = FALSE)

aa<-rownames(dev)
#write.table(as.matrix(aa),file="/data/kuw/biocore/wlku/Keji/wbc_285_read1_mapped/conet_het_result/conet_het_scrna_sch3k4_rhigh_het_clow_het2_motif_mono.txt",sep="\t",quote=FALSE,row.names=TRUE, na = "1000")
write.table(as.matrix(aa),file="./data/output/Figure5/conet_het_result/conet_het_scrna_sch3k4_rhigh_het_clow_het2_motif_mono.txt",sep="\t",quote=FALSE,row.names=TRUE, na = "1000")

aa<-assays(dev)$z
#write.table(data.frame(aa),file="/data/kuw/biocore/wlku/Keji/wbc_285_read1_mapped/chromvar/conet_het_scrna_sch3k4_rhigh_het_clow_het2_zscore_mono.txt",sep="\t",quote=FALSE,row.names=TRUE,na = "1000")
write.table(data.frame(aa),file="./data/output/Figure5/chromvar/conet_het_scrna_sch3k4_rhigh_het_clow_het2_zscore_mono.txt",sep="\t",quote=FALSE,row.names=TRUE,na = "1000")

aa<-assays(dev)$deviations
#write.table(data.frame(aa),file="/data/kuw/biocore/wlku/Keji/wbc_285_read1_mapped/chromvar/conet_het_scrna_sch3k4_rhigh_het_clow_het2_dev_mono.txt",sep="\t",quote=FALSE,row.names=TRUE,na = "1000")
#write.table(data.frame(variability),file="/data/kuw/biocore/wlku/Keji/wbc_285_read1_mapped/chromvar/conet_het_scrna_sch3k4_rhigh_het_clow_het2_var_mono.txt",sep="\t",quote=FALSE,row.names=TRUE,na = "1000")
write.table(data.frame(aa),file="./data/output/Figure5/chromvar/conet_het_scrna_sch3k4_rhigh_het_clow_het2_dev_mono.txt",sep="\t",quote=FALSE,row.names=TRUE,na = "1000")
write.table(data.frame(variability),file="./data/output/Figure5/chromvar/conet_het_scrna_sch3k4_rhigh_het_clow_het2_var_mono.txt",sep="\t",quote=FALSE,row.names=TRUE,na = "1000")

