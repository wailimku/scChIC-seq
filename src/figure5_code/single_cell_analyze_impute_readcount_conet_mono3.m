a= readtable('../../data/input/Figure3/bulk_wbc_rc_52798_72_mat.txt');
bulk_peakname= table2array(a(:,1));
bulk_rc = table2array(a(:,2:73));

bulk_cpm = bulk_rc;
for i = 1:72
    bulk_cpm(:,i) = log2(bulk_rc(:,i)*1000000/sum(bulk_rc(:,i))+1);
end;

a= readtable('../../data/input/Figure3/bulk_wbc_file_name.txt');
bulkfile = table2array(a(:,2));

wbc_cpm=bulk_cpm; 
wbc_rc = bulk_rc;

wbc_cpm_type = zeros(max(size(wbc_cpm)),12);
%wbc_cpm_type(:,1) = mean(wbc_cpm(:,1:10),2);
wbc_cpm_type(:,1) = mean(wbc_cpm(:,3:7),2);
wbc_cpm_type(:,2) = mean(wbc_cpm(:,11:16),2);
wbc_cpm_type(:,3) = mean(wbc_cpm(:,17:24),2);
wbc_cpm_type(:,4) = mean(wbc_cpm(:,25:30),2);
wbc_cpm_type(:,5) = mean(wbc_cpm(:,31:47),2);
wbc_cpm_type(:,6) = wbc_cpm(:,48);
wbc_cpm_type(:,7) = mean(wbc_cpm(:,49:50),2);
%wbc_cpm_type(:,7) = mean(wbc_cpm(:,[17:24,35:30,31:47,49,50,62]),2);

wbc_cpm_type(:,8) = mean(wbc_cpm(:,51:59),2);
wbc_cpm_type(:,9) = mean(wbc_cpm(:,60:61),2);
wbc_cpm_type(:,10) = wbc_cpm(:,62);
wbc_cpm_type(:,11) = mean(wbc_cpm(:,63:64),2);
wbc_cpm_type(:,12) = mean(wbc_cpm(:,65:72),2);

wbc_rc_type = zeros(max(size(wbc_rc)),12);
%wbc_rc_type(:,1) = mean(wbc_rc(:,1:10),2);
wbc_rc_type(:,1) = mean(wbc_rc(:,3:7),2);
wbc_rc_type(:,2) = mean(wbc_rc(:,11:16),2);
wbc_rc_type(:,3) = mean(wbc_rc(:,17:24),2);
wbc_rc_type(:,4) = mean(wbc_rc(:,25:30),2);
wbc_rc_type(:,5) = mean(wbc_rc(:,31:47),2);
wbc_rc_type(:,6) = wbc_rc(:,48);
wbc_rc_type(:,7) = mean(wbc_rc(:,49:50),2);
wbc_rc_type(:,8) = mean(wbc_rc(:,51:59),2);
wbc_rc_type(:,9) = mean(wbc_rc(:,60:61),2);
wbc_rc_type(:,10) = wbc_rc(:,62);
wbc_rc_type(:,11) = mean(wbc_rc(:,63:64),2);
wbc_rc_type(:,12) = mean(wbc_rc(:,65:72),2);

wbc_rc_type(:,[3,4,6,7,8,9,12])=[];

wbc_cpm_type(:,[3,4,6,7,8,9,12])=[];

%================================%================================

a= readtable('../../data/input/Figure4/scimpute_Thscrna/scimpute_count.txt');
Th_gene = table2array(a(:,1));
Th_cnt = table2array(a(:,2:191));
Th_exp = Th_cnt;
for i = 1:190
    Th_exp(:,i) = log2(Th_exp(:,i)*1000000/sum(Th_exp(:,i))+1);
end;
Th_het = std(Th_exp')./mean(Th_exp');


a= readtable('../../data/input/Figure4/scimpute_Bscrna/scimpute_count.txt');
%a= readtable('/data/kuw/biocore/wlku/Keji/scimpute/B_cell_seurat_norm_data.txt');
B_gene = table2array(a(:,1));
B_cnt = table2array(a(:,2:331));
B_exp = B_cnt;

for i = 1:330
    B_exp(:,i) = log2(B_exp(:,i)*1000000/sum(B_exp(:,i))+1);
end;
B_het = std(B_exp')./mean(B_exp');


a= readtable('../../data/input/Figure4/scimpute_NKscrna/scimpute_count.txt');
%a= readtable('/data/kuw/biocore/wlku/Keji/scimpute/NK_cell_seurat_norm_data.txt');
NK_gene = table2array(a(:,1));
NK_cnt = table2array(a(:,2:769));
NK_exp = NK_cnt;
for i = 1:768
    NK_exp(:,i) = log2(NK_exp(:,i)*1000000/sum(NK_exp(:,i))+1);
end;
NK_het = std(NK_exp')./mean(NK_exp');


a= readtable('../../data/input/Figure4/scimpute_Monoscrna/scimpute_count.txt');
%a= readtable('/data/kuw/biocore/wlku/Keji/scimpute/Mono_seurat_norm_data.txt');
Mono_gene = table2array(a(:,1));
Mono_cnt = table2array(a(:,2:130));
Mono_exp = Mono_cnt;
for i = 1:129
    Mono_exp(:,i) = log2(Mono_exp(:,i)*1000000/sum(Mono_exp(:,i))+1);
end;
Mono_het = std(Mono_exp')./mean(Mono_exp');
%================================%================================

a= readtable('../../data/input/Figure4/t-atac/GSE107817_RAW/primary_single/tatac-seq_cd4_c6_peak_id.txt','Delimiter','\t');
%a= readtable('/data/kuw/biocore/wlku/Keji/t-atac/GSE107817_RAW/primary_single/Tatac-seq_all_peak_id.txt','Delimiter','\t');
cd4_atac_peak = table2array(a(:,2));

cd4_atac_gene = cd4_atac_peak;
cd4_dis = zeros(max(size(cd4_atac_gene)),1);
for i = 1:max(size(cd4_atac_peak))
    a = strsplit(char(cd4_atac_peak(i)),'_');
    cd4_atac_gene(i) = cellstr(a(5));
    cd4_dis(i) = str2num(char(a(4)));
end;


qdel = find(abs(cd4_dis)>0);
cd4_dis(qdel)=[];
cd4_atac_peak(qdel)=[];
cd4_atac_gene(qdel)=[];

a= readtable('../../data/input/Figure4/t-atac/GSE107817_RAW/primary_single/sci_impute/scimpute_count.txt');
%a= readtable('/data/kuw/biocore/wlku/Keji/t-atac/GSE107817_RAW/primary_single/sci_impute_t_atac_all/scimpute_count.txt');
cd4_atac_cnt = table2array(a(:,2:97));
cd4_atac_cnt(qdel,:)=[];
cd4_atac_exp= cd4_atac_cnt;
for i = 1:96
    cd4_atac_exp(:,i) = log2(cd4_atac_cnt(:,i)*1000000/sum(cd4_atac_cnt(:,i))+1);
end;
q = find(sum(cd4_atac_cnt)<1000);
cd4_atac_exp(:,q)=[];
cd4_atac_cnt(:,q)=[];
q = find(sum(cd4_atac_cnt')==0);
cd4_atac_het = std(cd4_atac_exp')./mean(cd4_atac_exp');
%-------------------------------------------------------------------
%-------------------------------------------------------------------
%-------------------------------------------------------------------
a= readtable('../../data/input/Figure4/greenleaf_scatac_scRNA/supplementary_code/output/atac-seq_pbmc_mono_peak_id.txt','Delimiter','\t');
mono_atac_peak = table2array(a(:,2));

mono_atac_gene = mono_atac_peak;
mono_dis = zeros(max(size(mono_atac_gene)),1);
for i = 1:max(size(mono_atac_peak))
    a = strsplit(char(mono_atac_peak(i)),'_');
    mono_atac_gene(i) = cellstr(a(5));
    mono_dis(i) = str2num(char(a(4)));
end;

qdel = find(abs(mono_dis)>0);
mono_dis(qdel)=[];
mono_atac_peak(qdel)=[];
mono_atac_gene(qdel)=[];

a= readtable('../../data/input/Figure4/greenleaf_scatac_scRNA/supplementary_code/output/atac-seq_pbmc_mono_scImpute_count_mat2.txt');
mono_atac_cnt = table2array(a(:,2:65));
mono_atac_cnt(qdel,:)=[];
mono_atac_exp= mono_atac_cnt;
for i = 1:64
   mono_atac_exp(:,i) = log2(mono_atac_cnt(:,i)*1000000/sum(mono_atac_cnt(:,i))+1);
end;
q = find(sum(mono_atac_cnt)<1000);
mono_atac_exp(:,q)=[];
mono_atac_cnt(:,q)=[];
q = find(sum(mono_atac_cnt')==0);
mono_atac_het = std(mono_atac_exp')./mean(mono_atac_exp');



%========================================================================
%========================================================================
%========================================================================
a= readtable('../../data/input/Figure4/scimpute/scimpute_count.txt');
%a= readtable('/data/kuw/biocore/wlku/Keji/scwbc_rc_52798_242_mat.txt');
rc = table2array(a(:,2:243));
a= readtable('../../data/input/Figure4/scimpute/peakname2.txt','Delimiter','\t','ReadVariableNames',0);
peakname = table2array(a);
peakgene = peakname;
dis = zeros(max(size(peakgene )),1);
for i = 1:max(size(peakname))
    a = strsplit(char(peakname(i)),'_');
    peakgene(i) = cellstr(a(5));
    dis(i) = str2num(char(a(4)));
end;

%a = readtable('/data/kuw/biocore/wlku/Keji/scimpute/wbc_ident2.txt')
%clus = table2array(a);
rc(:,112)=[];
sc3_cellclus = readtable('../../data/input/Figure4/cell_clus_sc3_12779_wimpute_nobatch_7a.txt');


clus = table2array(sc3_cellclus(:,2));

cpm = rc;
for i = 1:size(rc,2)
    cpm(:,i) = rc(:,i)*1000000/sum(rc(:,i));
end;    
q = find(clus==7);
cpm33 = cpm(:,q);
rc33 = rc(:,q);
q1 = find(sum(cpm33')==0);
cpm33(q1,:)=[];
peakname2 = peakname;
wbc_cpm_type2 = wbc_cpm_type;
peakgene2 = peakgene;
peakname2(q1)=[];
wbc_cpm_type2(q1,:)=[];
peakgene2(q1)=[];
dis2 = dis;
dis2(q1)=[];
rc33(q1,:)=[];

q = find(abs(dis2)>0);
dis2(q)=[];
cpm33(q,:)=[];
peakgene2(q)=[];
rc33(q,:)=[];
peakname2(q)=[];
wbc_cpm_type2(q,:)=[];

%+==============================================================
%+==============================================================
[ia,ib,ic] = intersect(mono_atac_gene,Mono_gene);
scrna_mat = Mono_exp(ic,:);
scatac_mat = mono_atac_exp(ib,:);
mono_atac_peak2 = mono_atac_peak(ib);

[ia2,ib2,ic2] = intersect(ia, peakgene2,'stable');
scchip_mat = log2(cpm33(ic2,:)+1);
scrna_mat = scrna_mat(ib2,:);
scatac_mat =scatac_mat(ib2,:);
mono_atac_peak3 = mono_atac_peak2(ib2);
peakname3 = peakname2(ic2);

q = find(sum(scrna_mat')==0|sum(scatac_mat')==0|sum(scchip_mat')==0);
scrna_mat(q,:)=[];
scatac_mat(q,:)=[];
scchip_mat(q,:)=[];
peakname3(q)=[];



peakgene3 = ia2;
peakgene3(q)=[];
rhet = std(scrna_mat')./mean(scrna_mat');
ahet = std(scatac_mat')./mean(scatac_mat');
chet = std(scchip_mat')./mean(scchip_mat');

atac_net = squareform(pdist(scatac_mat,'correlation'),'tomatrix');
rna_net = squareform(pdist(scrna_mat,'correlation'),'tomatrix');
chip_net = squareform(pdist(scchip_mat,'correlation'),'tomatrix');


%*******************************
rr1 = corrcoef(log2(scrna_mat'+1));
rr2 = corrcoef(log2(scatac_mat'+1));
rr3 = corrcoef(log2(scchip_mat'+1));
rrzs1 = rr1;
rrzs2 = rr2;
rrzs3 = rr3;
zs1 = zscore(rr1);
zs2 = zscore(rr2);
zs3 = zscore(rr3);
q1 = find(zs1<=0);
q2 = find(zs2<=0);
q3 = find(zs3<=0);
zs1(q1)=0;
zs2(q2)=0;
zs3(q3)=0;
tic
for i = 1:max(size(rr3))-1 %5732
    for j = i+1:max(size(rr3)) % 5733
        rrzs1(i,j) = sqrt(zs1(i,j).^2+ zs1(j,i).^2);
        rrzs2(i,j) = sqrt(zs2(i,j).^2+ zs2(j,i).^2);
        rrzs3(i,j) = sqrt(zs3(i,j).^2+ zs3(j,i).^2);
        rrzs1(j,i) = rrzs1(i,j);
        rrzs2(j,i) = rrzs2(i,j);        
        rrzs3(j,i) = rrzs3(i,j);
    end;
    toc;
end;   
rr1 = rrzs1;
rr2 = rrzs2;
rr3 = rrzs3;

rrc1 = rr1;
q = find(rrc1<1.65);
rrc1(q)=0;
q = find(rrc1>=1.65);
rrc1(q)=1;   

rrc2 = rr2;
q = find(rrc2<1.65);
rrc2(q)=0;
q = find(rrc2>=1.65);
rrc2(q)=1;   


rrc3 = rr3;
q = find(rrc3<1.65);
rrc3(q)=0;
q = find(rrc3>=1.65);
rrc3(q)=1;   

G1 = graph(rrc1);
G2 = graph(rrc2);
G3 = graph(rrc3);

pg_rank1 = centrality(G1, 'Degree');
pg_rank2 = centrality(G2, 'Degree');
pg_rank3 = centrality(G3, 'Degree');

xx0=[pg_rank1,pg_rank2,pg_rank3,rhet',ahet',chet'];
xx=[pg_rank1,pg_rank2,pg_rank3,rhet',ahet',chet'];
%xx=[pg_rank1,pg_rank3,rhet',chet'];
%xx=[pg_rank1,pg_rank2,rhet',ahet'];
xx2 = quantilenorm(xx);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sel1 = 4;
sel2 = 1;
sel3 = 3;
sel4 = 6;
sel5 = 5;
sel6 = 2;

q = find(xx2(:,sel6)./xx2(:,sel5)<2 & xx2(:,sel5)./xx2(:,sel6)<2 & xx2(:,sel1)./xx2(:,sel2)<2 & xx2(:,sel2)./xx2(:,sel1)<2 & xx2(:,sel3)./xx2(:,sel4)>1.5);
%q = find(xx2(:,sel6)./xx2(:,sel5)<1.5 & xx2(:,sel5)./xx2(:,sel6)<1.5 & xx2(:,sel1)./xx2(:,sel2)<1.5 & xx2(:,sel2)./xx2(:,sel1)<1.5 & xx2(:,sel3)./xx2(:,sel4)>1.5);
%biva1 = [H3K27me3_cpm3(q),H3K4me3_cpm3(q)];
peakname4 = peakname3(q);
%%writetable(table(peakname4),'/data/kuw/biocore/wlku/Keji/wbc_285_read1_mapped/conet_het_result/conet_het_scrna_sch3k4_rhigh_het_clow_het_tcell2.txt');
%writetable(table(peakname4),'/data/kuw/biocore/wlku/Keji/wbc_285_read1_mapped/conet_het_result/conet_het_scrna_sch3k4_rhigh_het_clow_het2_mono.txt');
writetable(table(peakname4),'../../data/output/Figure5/conet_het_result/conet_het_scrna_sch3k4_rhigh_het_clow_het2_mono.txt');

q = find(xx2(:,sel6)./xx2(:,sel5)<2 & xx2(:,sel5)./xx2(:,sel6)<2 & xx2(:,sel1)./xx2(:,sel2)<2 & xx2(:,sel2)./xx2(:,sel1)<2 & xx2(:,sel4)./xx2(:,sel3)>1.5);
%q = find(xx2(:,sel6)./xx2(:,sel5)<1.5 & xx2(:,sel5)./xx2(:,sel6)<1.5 & xx2(:,sel1)./xx2(:,sel2)<1.5 & xx2(:,sel2)./xx2(:,sel1)<1.5 & xx2(:,sel4)./xx2(:,sel3)>1.5);
%biva2 = [H3K27me3_cpm3(q),H3K4me3_cpm3(q)];
peakname4 = peakname3(q);
%%writetable(table(peakname4),'/data/kuw/biocore/wlku/Keji/wbc_285_read1_mapped/conet_het_result/conet_het_scrna_sch3k4_rhigh_het_chigh_het_tcell2.txt');
%writetable(table(peakname4),'/data/kuw/biocore/wlku/Keji/wbc_285_read1_mapped/conet_het_result/conet_het_scrna_sch3k4_rhigh_het_chigh_het2_mono.txt');
writetable(table(peakname4),'../../data/output/Figure5/conet_het_result/conet_het_scrna_sch3k4_rhigh_het_chigh_het2_mono.txt');



