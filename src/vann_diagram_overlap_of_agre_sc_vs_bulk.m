

path1 = '../data/input/';
file_name1 = 'singlecell_read1-W200-G200-E.01_with_peak_pos_sorted.txt';

    
path2 = '../data/input/';
file_name2 = 'white_blood_cell_sicer_peaks_combined_peaks_sorted_combined.txt';

fp = fopen(strcat(char(path1), char(file_name1)), 'r');

nlines = 0;
while~feof(fp)
    a = fscanf(fp, '%s', 5);
    nlines = nlines + 1;
end;
fclose(fp);
nlines = nlines - 1;


chr1 =  cell([nlines, 1]);
index1 = zeros(nlines, 1);
peak_ss1 = zeros(nlines, 1);
peak_es1 = zeros(nlines, 1);
peak_pos1 = zeros(nlines, 1);


fp = fopen(strcat(char(path1), char(file_name1)), 'r');

for i =1: nlines
    a = fscanf(fp, '%s', 1);
    chr1(i) = cellstr(a);
    index1(i) = fscanf(fp, '%d', 1);
 
    b = fscanf(fp, '%d', 1);
    peak_ss1(i) = b;
    b = fscanf(fp, '%d', 1);
    peak_es1(i) = b;
    peak_pos1(i) = fscanf(fp, '%d', 1);
end;    
fclose(fp);

q= [];
for i=2:nlines-1
    if(peak_ss1(i)-peak_es1(i-1)<1000 | peak_ss1(i+1)-peak_ss1(i)<1000)
        q=[q,i];
    end;
end;    

chr1(q)=[];
index1(q)=[];
peak_ss1(q)=[];
peak_es1(q)=[];
peak_pos1(q)=[];

fp = fopen(strcat(char(path1),'singlecell_read1-W200-G200-E.01_filtered_1000_with_peak_pos_sorted.txt'),'w');
for i = 1:max(size(chr1))
    fprintf(fp,'%s\t%d\t%d\t%d\t%d\n',char(chr1(i)), index1(i),peak_ss1(i),peak_es1(i),peak_pos1(i));
end;
fclose(fp);

%======================================================================================
fp = fopen(strcat(char(path2), char(file_name2)), 'r');

nlines = 0;
while~feof(fp)
    a = fscanf(fp, '%s', 5);
    nlines = nlines + 1;
end;
fclose(fp);
nlines = nlines - 1;

chr2 =  cell([nlines, 1]);
peak_ss2 = zeros(nlines, 1);
peak_es2 = zeros(nlines, 1);

fp = fopen(strcat(char(path2), char(file_name2)), 'r');


for i =1: nlines
    a = fscanf(fp, '%s', 1);
    chr2(i) = cellstr(a);
    a = fscanf(fp, '%s', 1);
    b = fscanf(fp, '%d', 1);
    peak_ss2(i) = b;
    b = fscanf(fp, '%d', 1);
    peak_es2(i) = b;
    b = fscanf(fp, '%d', 1);
end;    
fclose(fp);


%=========================================================================

test_chr1 = chr1 ;
test_chr2 = chr2 ;
size_chr1 = max(size(chr1));
size_chr2 = max(size(chr2));
overlap_chr1 = zeros(size_chr1, 1);

tic;
for i = 1: size_chr1
    q = find(strcmp(chr2, chr1(i))==1);
    if(min(size(q))>0)
        for k = 1: max(size(q))
            j = q(k);
            if(peak_ss1(i)>=peak_ss2(j) & peak_ss1(i)<=peak_es2(j))
                overlap_chr1(i) = overlap_chr1(i) +1;
                test_chr1(i) = chr2(j);
            elseif(peak_es1(i)>=peak_ss2(j) & peak_es1(i)<=peak_es2(j))
                overlap_chr1(i) = overlap_chr1(i) +1;  
                test_chr1(i) = chr2(j);
            elseif(peak_ss2(j)>=peak_ss1(i) & peak_ss2(j)<=peak_es1(i))        
                overlap_chr1(i) = overlap_chr1(i) +1;
                test_chr1(i) = chr2(j);
            elseif(peak_es2(j)>=peak_ss1(i) & peak_es2(j)<=peak_es1(i))         
                overlap_chr1(i) = overlap_chr1(i) +1;
                test_chr1(i) = chr2(j);
            end;
        end;
    end;    
    clear q;
    %if(floor(i/2000)*2000==i)
        toc
    %end;    
end;   

overlap_chr2 = zeros(size_chr2, 1);

tic;
for i = 1: size_chr2
    q = find(strcmp(chr1, chr2(i))==1);
    if(min(size(q))>0)
        for k = 1: max(size(q))
            j = q(k);
            if(peak_ss2(i)>=peak_ss1(j) & peak_ss2(i)<=peak_es1(j))
                overlap_chr2(i) = overlap_chr2(i) +1;
                test_chr2(i) = chr1(j);
            elseif(peak_es2(i)>=peak_ss1(j) & peak_es2(i)<=peak_es1(j))
                overlap_chr2(i) = overlap_chr2(i) +1;
                test_chr2(i) = chr1(j);
            elseif(peak_ss1(j)>=peak_ss2(i) & peak_ss1(j)<=peak_es2(i))        
                overlap_chr2(i) = overlap_chr2(i) +1;
                test_chr2(i) = chr1(j);
            elseif(peak_es1(j)>=peak_ss2(i) & peak_es1(j)<=peak_es2(i))         
                overlap_chr2(i) = overlap_chr2(i) +1;
                test_chr2(i) = chr1(j);
            end;
        end;
    end;    
    clear q;

    %if(floor(i/2000)*2000==i)
        toc
    %end;       
end;    

%clear peak_ss1 peak_ss2 peak_es1 peak_es2;
overlap_chr = zeros(4,1);
%==========================================================================

for i = 1: size_chr2
    if(overlap_chr2(i)>0 & strcmp(chr2(i), test_chr2(i))==1)
        overlap_chr(3) = overlap_chr(3) + 1;    
    end;
end;    


for i = 1: size_chr2
    if(overlap_chr2(i)==0)
        overlap_chr(4) = overlap_chr(4) + 1;
    end;
end;    

for i = 1: size_chr1
    if(overlap_chr1(i)>0 & strcmp(chr1(i), test_chr1(i))==1)
        overlap_chr(1) = overlap_chr(1) + 1;     
    end;
end;    

for i = 1: size_chr1
    if(overlap_chr1(i)==0)
        overlap_chr(2) = overlap_chr(2) + 1;     
    end;
end;    

%clear test_chr1 test_chr2 overlap_chr1 overlap_chr2 chr1 chr2;
%clearvars -except overlap_chr;
fprintf('%s\t %s\t %s\t %s\n', 'overlap', 'non-overlap', 'total','fraction of overlap');
fprintf('%d\t %d\t %d\t %f\n', overlap_chr(1), overlap_chr(2), overlap_chr(1)+overlap_chr(2), overlap_chr(1)/(overlap_chr(1)+overlap_chr(2)));
fprintf('%d\t %d\t %d\t %f\n', overlap_chr(3), overlap_chr(4), overlap_chr(3)+overlap_chr(4), overlap_chr(3)/(overlap_chr(3)+overlap_chr(4)));



