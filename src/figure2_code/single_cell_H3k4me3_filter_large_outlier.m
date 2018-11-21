a = readtable('../../data/input/Figure2/wc_sc1_bed.txt','ReadvariableNames',0,'delimiter','\t');
b1 = table2array(a(:,1));

filename = table2array(a(:,2));
b3 = log2(b1);
[b4,ib] = sort(b3);
filename = filename(ib);
q = zeros(285,1);

for i = 1:10:280  
    q1 =isoutlier(b4,'movmedian',i);
    q = q + q1;
end;

q2 = find(q>0& b4>mean(b4));


filename(q2)=[];

fp = fopen('../../GSE105012/script_cat_sc','w');
fprintf(fp,'%s\t','cat')
for i = 1:max(size(filename))
    fprintf(fp,'%s\t',char(filename(i)));
end;    
fprintf(fp,'%s\t%s\n','>','combined_sc_white_blood_cell.bed');
fclose(fp);
