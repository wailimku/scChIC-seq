path = '../../GSE105012/';
listing = dir(strcat(char(path),'*noDup.bed.txt'));
tic
msize = max(size(listing));
filename = cell([msize, 1]);
for i = 1:msize
    filename(i) = cellstr(listing(i).name);
    toc
end;    

%====================================================================
%====================================================================
path2 = '../../data/temp/Figure2/';
peakfile = 'combined_sc_white_blood_cell-W200-G200-E.01.scoreisland';
fp = fopen(strcat(char(path2), char(peakfile)), 'r');

nlines = 0;
while~feof(fp)
    a = fscanf(fp, '%s', 4);
    nlines = nlines + 1;
end;
fclose(fp);
nlines = nlines - 1;
chr1 =  cell([nlines, 1]);
peak_ss1 = zeros(nlines, 1);
peak_es1 = zeros(nlines, 1);
fp = fopen(strcat(char(path2), char(peakfile)), 'r');

for i =1: nlines
    a = fscanf(fp, '%s', 1);
    chr1(i) = cellstr(a);
    b = fscanf(fp, '%d', 1);
    peak_ss1(i) = b + 2000;
    b = fscanf(fp, '%d', 1);
    peak_es1(i) = b - 2000;
    b = fscanf(fp, '%f', 1);
end;    
fclose(fp);
%====================================================================
%====================================================================

uni_chr = unique(chr1);
tic
for i =1:max(size(filename))
    a = strsplit(char(filename(i)),'.txt');
    fileout = strcat(a(1),'_filtered.txt');
    read1 = readtable(strcat(char(path), char(filename(i))));
    nsize = max(size(read1));
    chr2 = table2array(read1(:,1));
    read_ss = table2array(read1(:,2));
    read_es = table2array(read1(:,3));
    read_pos = floor((read_ss+read_es)/2);
    checkpos = zeros(nsize, 1);

    for j = 1:24
        q = find(strcmp(chr2, uni_chr(j)));
        q2 = find(strcmp(chr1, uni_chr(j)));
        for k = 1:max(size(q2))
            q3 = find(read_pos(q)<=peak_es1(q2(k)) & read_pos(q)>=peak_ss1(q2(k)));
            checkpos(q(q3))=1;
        end;    
    end;
    q4 = find(checkpos==1);
    qmis = setdiff(1:1:max(size(read_ss)), q4);
    read1(qmis,:)=[];
    writetable(read1,strcat('../../data/temp/Figure2/filtered_bed/',char(fileout)),'Delimiter','\t','WriteVariableNames',0);
end;    



