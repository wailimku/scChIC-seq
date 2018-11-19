# CRK (updating 18-Nov-2018)



A. Download files
--------------------------------------

1. Get the files from github and change directory to CRK

<pre>
[wlku@matrix] git clone https://github.com/wailimku/CRK.git
[wlku@matrix] cd CRK
[wlku@matrix CRK] ls
[wlku@matrix CRK] AdvancedColormap.m  <b>data</b>  <b>Figures</b>  <b>GSE105012</b>  README.md  <b>src</b> violin.m  violinplot.m
</pre>

2. Download the GSE105012.tar from GEO and put it into the folder  <b>GSE105012</b>
<pre>
[wlku@matrix CRK] cd GSE105012
[wlku@matrix GSE105012] ls
[wlku@matrix GSE105012] GSE105012.tar
[wlku@matrix GSE105012] tar -xvf GSE105012.tar
[wlku@matrix GSE105012] gunzip *.gz
[wlku@matrix GSE105012] cd ..
</pre>



B. H3K4me3 profiling  using low cell number 3T3 cells
--------------------------------------

1. Call peaks using SICER [1] for libraries using 3000, 1000, 300 and 100 cells. Peak files are savd in the directory './data/temp/Figure1/'
<pre>
[wlku@matrix CRK] sh ./src/figure1_code/script_sicer_low_cell
</pre>

2. Compute the TSS density profiles using HOMER [2] and plot the TSS profile plots using matlab (<b>Figure 1b</b>)
<pre>
[wlku@matrix CRK] sh ./src/figure1_code/script_run_homer_1
[wlku@matrix CRK] matlab -nodesktop
>> run ./src/Figure1_code/single_cell_3T3_tss_profile_plots(1,0)
</pre>

<p>A figure (Figure 1b) is shown by matlab</p>
<p><img src="Figures/Figure1/3T3_TSS_profiles.jpg" alt="Fig 1b" width="400" vspace="20" hspace="5" >
</p>

3. Plot the scatter plots between 3T3 bulk cell ChIP-data seq data and low cell number data
<pre>
[wlku@matrix CRK] Rscript ./src/Figure1_code/single_cell_analysze_H3K4me3_low_cell_num_scatter_bed2table.r
[wlku@matrix CRK] Rscript ./src/Figure1_code/single_cell_analysze_H3K4me3_low_cell_num_scatter_plot.r temp 1
[wlku@matrix CRK] cd Figures/Figure1/
[wlku@matrix CRK] ls
[wlku@matrix CRK] scatter_3T3_3000cell_at_H3K4me3_peaks_Fig1.pdf
[wlku@matrix CRK] xpdf scatter_3T3_3000cell_at_H3K4me3_peaks_Fig1.pdf
</pre>

<p><img src="Figures/Figure1/scatter_3T3_3000cell_at_H3K4me3_peaks_Fig1.jpg" alt="Fig 1d" width="500" vspace="20">
</p>
<div style="text-align:justify;">
This is Figure 1d
</div>

<p></p>

4. Compute the peak overlap between Bulk Cell ChIP-seq and low cell number library
<pre>
[wlku@matrix CRK] cd ../../
[wlku@matrix CRK] Rscript ./src/Figure1_code/peaks_overlap_for_3T3_cells.r 1
[1] "Bulk cell, number of peaks = "
[1] 16773
[1] "3000 cells scChIC-seq, number of peaks = "
[1] 16691
[1] "Overlap_peaks = "
[1] 13728
</pre>


5. Compute Clustering heatmap for 3T3, ESC, and Naive T cells
<pre>
[wlku@matrix CRK] matlab -nodesktop
>> run ./src/Figure1_code/analyze_H3K4me3_diff_cells_w_routput
</pre>
<p></p>

<p>A figure (Figure 1e) is shown by matlab</p>

<p></p> 
C. H3K4me3 profiling in Human white blood cells (WBCs)
--------------------------------------
1. Compute Clustering heatmap for 3T3, ESC, and Naive T cells



References:
--------------------------------------

1. SICER
2. HOMER
