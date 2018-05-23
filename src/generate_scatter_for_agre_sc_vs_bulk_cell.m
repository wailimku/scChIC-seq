%========================================
a = load('../data/output/single_cell_rc_data.mat');
ff = a.ff;
sc_peak_id = a.sc_peak_id;
chr = a.chr;
peakss = a.peakss;
peakes = a.peakes;
rc = a.rc;
cpm = a.cpm;
ib_spec  =a.ib_spec;
ic_spec = a.ic_spec;
clear a;
%========================================
a = load('../data/input/H3K4me3_wbc_ChIP-seq_ENCODE2.mat');
wbc_id = a.wbc_id;
wbc_rc = a.wbc_rc;
wbc_cpm_type = a.wbc_cpm_type;
wbc_rc_type = a.wbc_rc_type;
wbc_type_name = a.wbc_type_name;
clear a;
%========================================



cpm2 = sum(rc')*1000000/sum(sum(rc))';

wbc_cpm2 = sum(wbc_rc')*1000000/sum(sum(wbc_rc))';

ffig2 = figure;
scatter(log2(cpm2+1), log2(wbc_cpm2+1),'filled');
xlim([2,12]);
ylim([2,12]);
ylabel('Bulk cells');
xlabel('Pooled single cells');
cc = corrcoef(log2(cpm2+1),log2(wbc_cpm2+1));
title(strcat('P.C. = ', num2str(cc(1,2))));
saveas(ffig2,'../Figure/Figure2C_scatter','pdf');