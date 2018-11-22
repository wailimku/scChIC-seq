texx = 0;

%a= readtable('/data/kuw/biocore/wlku/Keji/bulk_wbc_rc_52798_72_mat.txt');
a= readtable('../../data/input/Figure3/bulk_wbc_rc_52798_72_mat.txt');

bulk_peakname= table2array(a(:,1));
bulk_rc = table2array(a(:,2:73));

bulk_cpm = bulk_rc;
for i = 1:72
    bulk_cpm(:,i) = bulk_rc(:,i)*1000000/sum(bulk_rc(:,i));
end;
a= readtable('../../data/input/Figure3/bulk_wbc_file_name.txt');
bulkfile = table2array(a(:,2));
wbc_type_name={'B_cell',
'CD14-positive_monocyte',
'CD4-positive_CD25-positive_alpha-beta_regulatory_T_cell',
'CD4-positive_alpha-beta_memory_T_cell',
'CD4-positive_helper_T_cell',
'CD8-positive_alpha-beta_T_cell',
'T-cell',
'common_myeloid_progenitor_CD34-positive',
'mononuclear_cell',
'naiveT',
'natural_killer_cell',
'neutrophil'};
wbc_cpm=log2(bulk_cpm+1);  
bulk_cpm = wbc_cpm;
%bulk_cpm = quantilenorm(bulk_cpm);
wbc_cpm = bulk_cpm;

wbc_rc = bulk_rc;
wbc_cpm_type = zeros(max(size(wbc_cpm)),12);
wbc_cpm_type(:,1) = mean(wbc_cpm(:,3:7),2);
wbc_cpm_type(:,2) = mean(wbc_cpm(:,11:16),2);
wbc_cpm_type(:,3) = mean(wbc_cpm(:,17:24),2);
wbc_cpm_type(:,4) = mean(wbc_cpm(:,25:30),2);
wbc_cpm_type(:,5) = mean(wbc_cpm(:,35:47),2);
wbc_cpm_type(:,6) = wbc_cpm(:,48);
wbc_cpm_type(:,7) = mean(wbc_cpm(:,49:50),2);
wbc_cpm_type(:,8) = mean(wbc_cpm(:,51:59),2);
wbc_cpm_type(:,9) = mean(wbc_cpm(:,60:61),2);
wbc_cpm_type(:,10) = wbc_cpm(:,62);
wbc_cpm_type(:,11) = mean(wbc_cpm(:,61:62),2);
wbc_cpm_type(:,12) = mean(wbc_cpm(:,63:64),2);

wbc_cpm_type(:,[3,4,6,7,8,9])=[];

import bioma.data.DataMatrix     
DM_B = DataMatrix(bulk_cpm(:,3:7));%3-7
DM_Mono = DataMatrix(bulk_cpm(:,11:16));%11-16
DM_T = DataMatrix(bulk_cpm(:,35:47));% 31-47
DM_NK = DataMatrix(bulk_cpm(:,61:62));

[PValues_T_B, TScores1] = mattest(DM_T, DM_B);
[PValues_T_NK, TScores1] = mattest(DM_T, DM_NK);
[PValues_T_Mono, TScores1] = mattest(DM_T, DM_Mono);

[PValues_B_T, TScores1] = mattest(DM_B, DM_T);
[PValues_B_NK, TScores1] = mattest(DM_B, DM_NK);
[PValues_B_Mono, TScores1] = mattest(DM_B, DM_Mono);

[PValues_NK_B, TScores1] = mattest(DM_NK, DM_B);
[PValues_NK_T, TScores1] = mattest(DM_NK, DM_T);
[PValues_NK_Mono, TScores1] = mattest(DM_NK, DM_Mono);

[PValues_Mono_B, TScores1] = mattest(DM_Mono, DM_B);
[PValues_Mono_T, TScores1] = mattest(DM_Mono, DM_T);
[PValues_Mono_NK, TScores1] = mattest(DM_Mono, DM_NK);
FValues_T_B = mafdr(PValues_T_B);
FValues_T_NK = mafdr(PValues_T_NK);
FValues_T_Mono = mafdr(PValues_T_Mono);
FValues_B_T = mafdr(PValues_B_T);
FValues_B_NK = mafdr(PValues_B_NK);
FValues_B_Mono = mafdr(PValues_B_Mono);
FValues_NK_B = mafdr(PValues_NK_B);
FValues_NK_T = mafdr(PValues_NK_T);
FValues_NK_Mono = mafdr(PValues_NK_Mono);
FValues_Mono_B = mafdr(PValues_Mono_B);
FValues_Mono_T = mafdr(PValues_Mono_T);
FValues_Mono_NK = mafdr(PValues_Mono_NK);

bulk_peakgene = bulk_peakname;
bulk_dis = zeros(max(size(bulk_peakname)),1);
for i = 1:max(size(bulk_peakname))
    a = strsplit(char(bulk_peakname(i)),'_');
    bulk_peakgene(i) = cellstr(a(5));
    bulk_dis(i) = str2num(char(a(4)));
end;

pcut=0.05;
fcut = 1.58;
ncut = log2(1+1);
n2cut = log2(3+1);
q_B = find(FValues_B_T<pcut & FValues_B_Mono<pcut & FValues_B_NK<pcut & wbc_cpm_type(:,1)-wbc_cpm_type(:,2)>fcut & wbc_cpm_type(:,1)-wbc_cpm_type(:,3)>fcut & wbc_cpm_type(:,1)-wbc_cpm_type(:,5)>fcut & wbc_cpm_type(:,1)>ncut &  wbc_cpm_type(:,2)<n2cut &  wbc_cpm_type(:,3)<n2cut &  wbc_cpm_type(:,5)<n2cut & abs(bulk_dis)<1);
q_Mono = find(FValues_Mono_T<pcut & FValues_Mono_B<pcut & FValues_Mono_NK<pcut & wbc_cpm_type(:,2)-wbc_cpm_type(:,1)>fcut & wbc_cpm_type(:,2)-wbc_cpm_type(:,3)>fcut & wbc_cpm_type(:,2)-wbc_cpm_type(:,5)>fcut & wbc_cpm_type(:,2)>ncut &  wbc_cpm_type(:,1)<n2cut &  wbc_cpm_type(:,3)<n2cut &  wbc_cpm_type(:,5)<n2cut  & abs(bulk_dis)<1);
q_T = find(FValues_T_B<pcut & FValues_T_Mono<pcut & FValues_T_NK<pcut & wbc_cpm_type(:,3)-wbc_cpm_type(:,1)>fcut & wbc_cpm_type(:,3)-wbc_cpm_type(:,2)>fcut & wbc_cpm_type(:,3)-wbc_cpm_type(:,5)>fcut & wbc_cpm_type(:,3)>ncut &  wbc_cpm_type(:,2)<n2cut &  wbc_cpm_type(:,1)<n2cut &  wbc_cpm_type(:,5)<n2cut  & abs(bulk_dis)<1);
q_NK = find(FValues_NK_T<pcut & FValues_NK_Mono<pcut & FValues_NK_B<pcut & wbc_cpm_type(:,5)-wbc_cpm_type(:,2)>fcut & wbc_cpm_type(:,5)-wbc_cpm_type(:,3)>fcut & wbc_cpm_type(:,5)-wbc_cpm_type(:,1)>fcut & wbc_cpm_type(:,5)>ncut &  wbc_cpm_type(:,2)<n2cut &  wbc_cpm_type(:,3)<n2cut &  wbc_cpm_type(:,1)<n2cut  & abs(bulk_dis)<1);

%===================================================================
sc3_para = '7a';
%sc3_result = readtable('/data/kuw/biocore/wlku/Keji/wbc_285_read1_mapped/SC3/clus_sc3_12779_wimpute_nobatch_7a.txt'); %4524-7c works
sc3_result = readtable('../../data/input/Figure3/clus_sc3_12779_wimpute_nobatch_7a.txt'); %4524-7c works
%sc3_result = readtable(strcat('./data/output/Figure3/SC3/clus_sc3_wimpute_nobatch_',char(sc3_para),'.txt')); %4524-7c works
sc3_clus_mark = cell2mat(table2cell(sc3_result));

%sc3_pval_result = readtable('./data/input/Figure3/pval_sc3_12779_wimpute_nobatch_7a.txt');
%sc3_pval_result = readtable('./data/output/Figure3/SC3/pval_sc3_wimpute_nobatch_7a.txt');

%sc3_pval_result = readtable('/data/kuw/biocore/wlku/Keji/wbc_285_read1_mapped/SC3/pval_sc3_12779_wimpute_nobatch_7a.txt');
sc3_pval_result = readtable('../../data/input/Figure3/pval_sc3_12779_wimpute_nobatch_7a.txt');

%sc3_pval_result = readtable(strcat('./data/output/Figure3/SC3/pval_sc3_wimpute_nobatch_',char(sc3_para),'.txt'));

sc3_pval = table2array(sc3_pval_result);
sc3_fdr = mafdr(sc3_pval(:,2));
%sc3_peak = readtable('/data/kuw/biocore/wlku/Keji/wbc_285_read1_mapped/SC3/peak_name_12779_wimpute_nobatch.txt','Delimiter','\t');
sc3_peak = readtable('../../data/input/Figure3/peak_name_12779_wimpute_nobatch.txt','Delimiter','\t');
%sc3_peak = readtable('./data/output/Figure3/SC3/peak_name_wimpute_nobatch.txt','Delimiter','\t');

sc3_peak_name = table2array(sc3_peak); 
[ia, qspec1, qspec] = intersect(bulk_peakname, sc3_peak_name,'stable');

%sc3_cellclus = readtable('/data/kuw/biocore/wlku/Keji/wbc_285_read1_mapped/SC3/cell_clus_sc3_12779_wimpute_nobatch_7a.txt');
sc3_cellclus = readtable('../../data/input/Figure3/cell_clus_sc3_12779_wimpute_nobatch_7a.txt');
%sc3_cellclus = readtable(strcat('./data/output/Figure3/SC3/cell_clus_sc3_wimpute_nobatch_',char(sc3_para),'.txt'));


sc3_cell = table2array(sc3_cellclus(:,2));

qp1 = find(sc3_pval(:,2)<=0.01);
%qp1 = find(sc3_fdr<=0.05);

uni_clus =unique(sc3_clus_mark(qp1,2));
clus_test = zeros(max(size(uni_clus)), 5);

qp1 = find(sc3_pval(:,2)<=1000);

tic
sco = zeros(max(size(uni_clus)),4);
for i = 1:max(size(uni_clus))
    qs1 = find(sc3_clus_mark(:,2)==uni_clus(i));
    
    %[qx,qy]=sort(sc3_pval(qs1,2));%0.00005%0.000000000 
    %[ia,ib] = sort(-log(sc3_pval(qs1,2)),'descend');
    %qx = ib(101:max(size(ia)));
    %qx = find(sc3_pval(qs1,2)>0.0000000000000005);
    %qx = find(sc3_pval(qs1,2)>0.00000000000000001);
    qx = find(sc3_fdr(qs1)>0.05);
    %qx = find(sc3_pval(qs1,2)>0.0000001);
    %qx = find(sc3_pval(qs1,2)>0.001);
    %qx = find(sc3_pval(qs1,2)>0.0015);
    %qx = find(sc3_fdr(qs1)>0.015);
    qs1(qx)=[];
    qnow = qspec1(qs1);
    %====
    %[ia2,ib2,ic] = intersect(bulk_peakgene(qnow),th_gene);
    %zz= -log2(sc3_pval(qs1(ib2),2));
    %rpkmnormdata(ic,:)
    %c1 = 7;
    %c2 = 2;
    %q1 = find(rpkmnormdata(ic,1)>c1 &  rpkmnormdata(ic,2)<c2 & rpkmnormdata(ic,3)<c2 & rpkmnormdata(ic,4)<c2)
    %q2 = find(rpkmnormdata(ic,2)>c1 &  rpkmnormdata(ic,1)<c2 & rpkmnormdata(ic,3)<c2 & rpkmnormdata(ic,4)<c2)
    %q3 = find(rpkmnormdata(ic,3)>c1 &  rpkmnormdata(ic,2)<c2 & rpkmnormdata(ic,1)<c2 & rpkmnormdata(ic,4)<c2)
    %q4 = find(rpkmnormdata(ic,4)>c1 &  rpkmnormdata(ic,2)<c2 & rpkmnormdata(ic,3)<c2 & rpkmnormdata(ic,1)<c2)
    
    %rpkmnormdata(ic(q1),:)
    %====
    %qnow = qspec1(qs1(qy(1:min(max(size(qy)),300))));
    
    q1=intersect(qnow,q_B);
    q2=intersect(qnow,q_Mono);
    q3=intersect(qnow,q_T);
    q4=intersect(qnow,q_NK);
    if(min(size(q1))>0)
        p1=1-sum(hygepdf(0:max(size(q1))*min(size(q1)),max(size(sc3_peak_name)),max(size((intersect(q_B,qspec1(qp1))'))),max(size(qnow))));
    else
        p1 = 1;
    end;
    if(min(size(q2))>0)
        p2=1-sum(hygepdf(0:max(size(q2))*min(size(q2)),max(size(sc3_peak_name)),max(size((intersect(q_Mono',qspec1(qp1))'))),max(size(qnow))));
    else
        p2 = 1;
    end;    
    if(min(size(q3))>0)    
        p3=1-sum(hygepdf(0:max(size(q3))*min(size(q3)),max(size(sc3_peak_name)),max(size((intersect(q_T,qspec1(qp1))'))),max(size(qnow))));
    else
        p3 = 1;
    end    
    if(min(size(q4))>0)
        p4=1-sum(hygepdf(0:max(size(q4))*min(size(q4)),max(size(sc3_peak_name)),max(size((intersect(q_NK,qspec1(qp1))'))),max(size(qnow))));
    else
        p4 = 1;
    end;    
    sco(i,1) = p1;
    sco(i,2) = p2;
    sco(i,3) = p3;
    sco(i,4) = p4;
end;    

fig1 = figure;
imagesc(-log2(sco));
set(gca,'YTick',[1,2,3,4,5],'YTicklabel',{'Cluster 1','Cluster 3','Cluster 4' ,'Cluster 6','Cluster 7'})
set(gca,'XTick',[1,2,3,4],'XTicklabel',{'B cell','Mono','T cell','NK cell'});
colorbar
caxis([0,6]);
%=========================================================================

t_marker_1 = 'chr6:30149600-30152599(RNF39)' ;
t_marker_2 = 'chr15:86599700-86602699(NTRK3)';
t_marker_3 = 'chr16:55179600-55182599(MT3)'  ;
nk_marker = 'chrX:50228700-50231699(DGKK)';
mono_marker1_1 ='chr5:38880100-38883099(OSMR)';
b_marker= 'chr1:239585600-239588599(RGS7)';
mono_marker2_1 ='chr12:50590400-50593399 (ACVRL1)';
mono_marker2_2 ='chr14:22356900-22359899(SLC7A7)';
mono_marker2_3 ='chr17:76896100-76899099(LINC00482)';



a = readtable('../../data/input/Figure3/consen_mat_sc3_12779_wimpute_nobatch_7a.txt');
%a = readtable('/data/kuw/biocore/wlku/Keji/wbc_285_read1_mapped/SC3/consen_mat_sc3_12779_wimpute_nobatch_7a.txt');
bb = table2array(a);
%sc3_cellclus = readtable('/data/kuw/biocore/wlku/Keji/wbc_285_read1_mapped/SC3/cell_clus_sc3_12779_wimpute_nobatch_7a.txt');
sc3_cellclus = readtable('../../data/input/Figure3/cell_clus_sc3_12779_wimpute_nobatch_7a.txt');

sc3_cell = table2array(sc3_cellclus(:,2));

%[v1,v2] = tsne(bb,'Algorithm','exact','NumPCAComponents',40,'Perplexity',200,'Distance','correlation');
[v1,v2] = tsne(bb,'Algorithm','exact','NumDimensions',3,'NumPCAComponents',144,'Perplexity',80,'Distance','correlation');



   
fig2 = Figure;
for i = 1:7
    q = find(sc3_cell==i);
    scatter3(v1(q,1),v1(q,2),v1(q,3),'filled');
    hold on
    if(texx==1)
        if(i==1)
            text(mean(v1(q,1)),mean(v1(q,2)), char(t_marker_1));
            text(mean(v1(q,1)),mean(v1(q,2)), char(t_marker_2));
            text(mean(v1(q,1)),mean(v1(q,2)), char(t_marker_3));
        end;  
        if(i==3)
            text(mean(v1(q,1)),mean(v1(q,2)), char(nk_marker));
        end;
        if(i==4)
            text(mean(v1(q,1)),mean(v1(q,2)), char(mono_marker1_1));
        end;
        if(i==6)
            text(mean(v1(q,1)),mean(v1(q,2)), char(b_marker));
        end;
        if(i==7)
            text(mean(v1(q,1)),mean(v1(q,2)), char(mono_marker2_1));
            text(mean(v1(q,1)),mean(v1(q,2)), char(mono_marker2_2));
            text(mean(v1(q,1)),mean(v1(q,2)), char(mono_marker2_3));
        end;    
    end;    
end; 

xlabel('tsne-1');
ylabel('tsne-2');
zlabel('tsne-3');

legend('T cell','No markers','NK cell', 'Monocytes','No markers', 'B cell', 'Monocytes');