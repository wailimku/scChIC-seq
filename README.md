# CRK


Update all by November-18-2018


-------------------------------------------------------------------

 last update: November-16-2018

 For Figure 1, Extended data Figs. 2-4
 In these figures, we generated H3K4me3 scChIC-seq data using low 
 cell number for NIH-3T3 cells, mESC cells, and Naïve T cells. We
 want to examine the specificity of profiling H3K4me3, and its 
 applicability to different cell types. 

-------------------------------------------------------------------


1) (This step is optional since you can directly download the output files from GEO data sets) Download the raw fastq data from Pubmed to the folder ./data/input/Figure1/. The list of files can be found in ./src/figure1_code/files_download.txt
	

2) (This step is optional since you can directly download the output files from PubMed) Run the script in /src/figure1_code/script_for_mapping for mapping. In the script file, the mapping command is 

a)	bowtie2 -p 18 -q -5 0 -3 0 -x /***/mm9/genome -U *.fastq |samtools view -bS - > *.bam	

Two output files, includes a bed file that contains unique reads and a bedgraph file for visualization on the genome browser, are for the downstream analyses.

3) Upload the bedgraph files to the washU genome browser[1].	 

4) (This step is optional since the files are prepared and saved in the folder ./src/input/Figure1/) Identify peaks using SICER

a)	sh ./src/Figure1_code/script_sicer
		

5) Plot TSS profiles plots (Figure 1b)

a)	Run script_run_homer_1 for calculating the tss density profile
b)	Run ./src/Figure1_code/single_cell_3T3_tss_profile_plots.m	

6) Plot scatter plots (Figure 1c and Extended Figure 4)

a)	Rscript ./src/Figure1_code/ single_cell_analysze_H3K4me3_low_cell_num_scatter_bed2table.r

b)	Rscript ./src/Figure1_code/single_cell_analysze_H3K4me3_low_cell_num_scatter_plot.r (Please change the parameter kkk for generating different figures)



7) Compute Overlap between peaks (Figure 1d and and Extended Figure 3c and 3d)
              
a)	Rscript ./src/Figure1_code/peaks_overlap_for_3T3_cells.r (Please change the parameter kkk for generating different figures)

8) Generate heatmap for 3T3, ESC, and Naïve T cells (Figure 1e)

a)	Run ./src/Figure1_code/ single_cell_analysze_H3K4me3_3T3_ESC_T_bed2table.r
b)	Run ./scr/Figure1_code/ analyze_log_log_plot_H3K4me3_diff_cells_w_routput.m

9) Generate TSS profiles plot for 3T3, ESC, and Naïve T cells	(Extended Data Fig. 4)

a) Run Homer to calculate the TSS profile
b) Run ./src/Figure1_code/single_cell_3T3_ESC_naive_tss_profile_plots.m

