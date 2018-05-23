%========================================
a = readtable('../data/input/tss_homer_output.txt');
tss = zeros(285,21);
for i = 1:285
    b = table2array(a(:,3+(i-1)*3));
    tss(i,:) = b';
end;
mtss = mean(tss)';
mm = max(mtss)/min(mtss);
index_tss = zeros(285,2);

for i = 1:285
    b= corrcoef(tss(i,:),mtss);
    index_tss(i,1) = b(1,2);
    if(min(tss(i,:))==0)
        index_tss(i,2) =max(tss(i,:)) - min(tss(i,:));
    else    
        index_tss(i,2) = (max(tss(i,:))/(min(tss(i,:))))/(mm*0.7);
    end;    
end;
qmis = find(index_tss(:,1)>0.9 & index_tss(:,2)>1);
qmis = setdiff(1:1:285,qmis);    
%========================================
a = load('../data/input/scfile_size_285.mat');
ff = a.ff;
mline = a.mline;
wbc_id = a.wbc_id_old;
path2 = '../../data/input/';
path2 = '../../wbc_285_read1_mapped/';
listing = dir(strcat(char(path2),'*_readcount_per_filtered_wbc.txt'));

msize = max(size(listing));
for i = 1:msize
    files1(i) = cellstr(listing(i).name);   
end;

qmis2 = find(mline<41427);
qmis3 = find(ff<4000 | ff>100000);

q = union(qmis3,union(qmis,qmis2));

files1(q)=[];
ff(q)=[];

msize = max(size(files1));
rc = zeros(41427,msize);
chr = cell([41427,msize]);
peakss = zeros(41427,msize);
peakes = zeros(41427,msize);

tic
for i = 1:msize
    b = readtable(strcat(char(path2), char(files1(i))));
    rc(:,i) = table2array(b(:,5));
    toc
end;
chr = table2array(b(:,1));
peakss = table2array(b(:,3));
peakes = table2array(b(:,4));


sc_peak_id = cell([41427,1]);
for i = 1:41427
    sc_peak_id(i) = cellstr(strcat(char(chr(i)),'_',num2str(peakss(i)),'_',num2str(peakes(i))));
end;    

[ia_spec,ib_spec,ic_spec] = intersect(sc_peak_id,wbc_id,'stable');

%+++++++++++++++++++++++++
%+++++++++++++++++++++++++
rc = rc(ib_spec,:);
chr = chr(ib_spec);
peakss = peakss(ib_spec);
peakes = peakes(ib_spec);
sc_peak_id = ia_spec;
wbc_id(ic_spec);
%+++++++++++++++++++++++++
%+++++++++++++++++++++++++
peak_length = peakes-peakss;
cpm = rc;
rc1 = rc;
wbc_node =1:1:max(size(cpm)); 
for i= 1:min(size(rc))
    q = find(rc(:,i)>=2);
    corr = cpm(:,i)*(ff(i)-sum(rc(:,i)))/20000;
    cpm(:,i) = (rc(:,i))*1000000*1000./((sum(rc(q,i)))*peak_length);
    rc1(:,i) = rc(:,i)*1000./peak_length;
end;

save('../data/output/single_cell_rc_data2','ff','sc_peak_id','chr','peakss','peakes','rc','cpm','ib_spec','ic_spec','files1');

