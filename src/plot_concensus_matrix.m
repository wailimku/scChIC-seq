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


%=============================
[qx,qy] = find(wbc_cpm_type>1);
qspec1 = unique(qx);
q = find(sum(cpm(qspec1,:)')==0 | std(log2(cpm(qspec1,:)+1)')./mean(log2(cpm(qspec1,:)+1)')<1);
qspec1(q)=[];

opt = statset('MaxIter',5000,'TolFun',0,'TolX',exp(-400),'Display','final');
[W,H] = nnmf(log2(cpm(qspec1,:)+1),25,'options',opt,'replicates',5);


batch_name = files1;
for i = 1:max(size(files1))
    s = strsplit(char(files1(i)),'_');
    batch_name(i) = cellstr(s(1));
end;    
uni_batch = unique(batch_name);
new_batch_vec = zeros(max(size(files1)), 1);

for i = 1:max(size(uni_batch))
    q = find(strcmp(uni_batch(i),batch_name)==1);
    new_batch_vec(q) = i;
end;    

batch_matrix = zeros(max(size(H)), 10);
for i =1:10
    ss = new_batch_vec;
    ss(:)=0.1;
    q = find(new_batch_vec==i);
    ss(q)=1;
    batch_matrix(:,i) = ss;
end;

test_H = zeros(25,10);
for i =1:25
    for j = 1:10
        a = corrcoef(batch_matrix(:,j),H(i,:)');
        test_H(i,j) = a(1,2);
    end;    
end;    
[qx,qy] = find(test_H>0.6);
qrem = setdiff(1:1:25,unique(qx));
rcpm6 = W(:,qrem)*H(qrem,:);
q = find(sum(rcpm6')>0);
rcpm7 = rcpm6(q,:);
%=========================================

cpm2 = rcpm7(:,sc0);
%cpm2 = cpm(:,sc0);
aa = squareform(pdist(cpm2','euclidean'),'tomatrix');
bb = exp(-aa/max(max(aa)));
deg = 1./(sum(bb).^0.5);

L =diag(ones(112,1))-diag(deg)*bb*diag(deg);
[v,w]=eigs(L,60);
concen_mat= zeros(112,112);

for i = 1:1000
    [ia,ib] = kmeans(corrcoef(v'),4);
    for j = 1:4
        q = find(ia==j);
        concen_mat(q,q) =  concen_mat(q,q) + 1;
    end;    
end;
concen_mat = concen_mat/1000;



%=========================================

a = readtable('consen_mat_sc3_8560.txt','delimiter',' ');