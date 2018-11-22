
a = readtable('../../data/temp/Figure2/output2.txt');

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
path2 = '../../data/temp/Figure2/filtered_bed/';
listing2 = dir(strcat(char(path2),'*filtered.txt'));


msize = max(size(listing2));
files2 = cell([msize, 1]);

for i = 1:msize
    files2(i) = cellstr(listing2(i).name);    
end;

tic
ff = zeros(msize, 1);
for i = 1:msize
    nline = 0;
    fp = fopen(strcat(char(path2), char(files2(i))),'r');
    while~feof(fp)   
        fgetl(fp);
        nline = nline + 1;
    end;
    fclose(fp);
    ff(i) = nline;
    toc
end;

qmis3 = find(ff<4000 | ff>100000);
q = union(qmis3,qmis);


files2(q)=[];
ff(q)=[];
t = table(files2);
t2 = table(ff);
writetable(t,'../../data/temp/Figure2/sel_242_file.txt','WriteVariableNames',0)
writetable(t2,'../../data/temp/Figure2/ff.txt','WriteVariableNames',0)

%==================================================================