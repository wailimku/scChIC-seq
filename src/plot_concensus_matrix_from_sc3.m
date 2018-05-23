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
a = load('../data/output/single_cell_rc_data2.mat');
files1 = a.files1;
clear a;
%========================================
a = load('../data/input/H3K4me3_wbc_ChIP-seq_ENCODE.mat');
wbc_id = a.wbc_id;
wbc_cpm_type = a.wbc_cpm_type;
wbc_rc_type = a.wbc_rc_type;
wbc_type_name = a.wbc_type_name;
clear a;
%========================================


sc3_cellclus = readtable('../data/input/cell_clus_sc3_8560_7a.txt');
sc3_cell = table2array(sc3_cellclus);
sc1 = find(sc3_cell(:,2)==6);
sc2 = find(sc3_cell(:,2)==7);
sc3 = find(sc3_cell(:,2)==4);
sc4 = find(sc3_cell(:,2)==5);

sc0 = [sc1',sc2',sc3',sc4'];

%=========================================
a = zeros(242, 242);
fp = fopen('../data/input/consen_mat_sc3_8560.txt','r');
b = fscanf(fp,'%s',242); 

for i = 1:242
    a(i,:) = fscanf(fp,'%f',242);
end;
imagesc(a(sc0,sc0));
AdvancedColormap('bwr');
caxis([0.4,0.8]);
colorbar