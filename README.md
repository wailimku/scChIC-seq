# CRK (updating 18-Nov-2018)



Download files
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



H3K4me3 profiles in 3T3 cells using low cell numer
--------------------------------------

1. Call peaks using SICER [1] for libraries using 3000, 1000, 300 and 100 cells. Peak files are savd in the directory './data/temp/Figure1/'
<pre>
[wlku@matrix CRK] sh ./src/figure1_code/script_sicer_low_cell
</pre>

2. Compute the TSS density profiles using HOMER [2] and plot the TSS profile plots using matlab (<b>Figure 1b</b>)
<pre>
[wlku@matrix CRK] sh ./src/figure1_code/script_run_homer_1
[wlku@matrix CRK] matlab -nodesktop
[wlku@matrix CRK] run ./src/Figure1_code/single_cell_3T3_tss_profile_plots(1,0)
</pre>

3. Plot the scatter plots between 3T3 bulk cell ChIP-data seq data and low cell number data
<pre>
[wlku@matrix CRK] Rscript ./src/Figure1_code/ single_cell_analysze_H3K4me3_low_cell_num_scatter_bed2table.r
</pre>

References:
--------------------------------------

1. SICER
2. HOMER
