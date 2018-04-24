

path1 = '/data/kuw/biocore/wlku/Keji/encode_white_blood_cell/';
listing = dir(strcat(char(path1), '*wbc.txt'));
wbc_chip_name = cell([max(size(listing)), 1]);

num_wbc_chip = max(size(listing));
for i = 1:num_wbc_chip
    wbc_chip_name(i) = cellstr(listing(i).name);
end;

nline = 0;
fp = fopen(strcat(char(path1),char(wbc_chip_name(1))),'r');
while~feof(fp)
    fgetl(fp);
    nline = nline + 1;
end;
fclose(fp);
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

wbc_rc = zeros(nline, num_wbc_chip);
wbc_chr = cell([num_wbc_chip, 1]);
wbc_id  = cell([num_wbc_chip, 1]);
wbc_peakss = zeros(num_wbc_chip, 1);
wbc_peakes = zeros(num_wbc_chip, 1);
tic;
for i=1:num_wbc_chip
    fp = fopen(strcat(char(path1),char(wbc_chip_name(i))), 'r');
    for j =1:nline
        a = fscanf(fp, '%s', 1);
        b1 = fscanf(fp, '%d', 1);
        b2 = fscanf(fp, '%d', 1);
        b3 = fscanf(fp, '%d', 1);
        wbc_rc(j,i) = fscanf(fp, '%d', 1);
        if(i==1)
            wbc_chr(j) = cellstr(a);
            wbc_peakss(j) = b2;
            wbc_peakes(j) = b3;
            wbc_id(j) = cellstr(strcat(char(a),'_',num2str(b2),'_',num2str(b3)));
        end;    
    end;
    fclose(fp);
    toc
end;

wbc_id = cell([41427,1]);
for i = 1:41427
    wbc_id(i) = cellstr(strcat(char(wbc_chr(i)),'_',num2str(wbc_peakss(i)),'_',num2str(wbc_peakes(i))));
end;    

peak_len = wbc_peakes - wbc_peakss;
wbc_cpm = wbc_rc;
wbc_rpkm = wbc_rc;
wbc_rc2=wbc_rc;
  
wbc_cpm = wbc_rc;
for i= 1:num_wbc_chip
    wbc_cpm(:,i) = wbc_rc(:,i)*100000*1000./(sum(wbc_rc(:,i))*peak_len);
end;
