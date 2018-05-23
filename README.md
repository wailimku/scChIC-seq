# CRK

Figure 1.

1. Go to src folder. run prepare_H3K4me3_bulk_cell_data.m

2. For TSS profile plot.
run single_cell_TSS_plot.m

3. Genrating the overlap between pooled single cells and bulk cells.
run vann_diagram_overlap_of_agre_sc_vs_bulk.m

4. (Optional, the file has been generated already in the output folder).
Otherwise, run filtering_reads_step_for_getting_242_sc.m

5. (Generate scatter plot), run generate_scatter_for_agre_sc_vs_bulk_cell.m

6. (Generate plots of precision and sensitivity)
run generate_pre_sen.m

Figure 2.

7. We computed the CPM for the 242 single cell H3K4me3 cut&run data. Remove the batch effect by using the non-negative matrix factorization. The final output of the corrected CPM and raw read count data are obtained as the input for SC3.

8. The consensus matrix are the output of SC3 (/data/input/consen_mat_sc3_8560.txt), the consensus matrix are plotted using the code "plot_concensus_matrix_from_SC3".

9.
