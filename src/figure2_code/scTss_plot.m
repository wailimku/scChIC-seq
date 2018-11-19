%=========================================
path1 = '../../data/temp/Figure2/';

a = readtable(strcat(char(path1),'tss_homer_output2.txt'));
tss = zeros(285,501);
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

fig1 = figure;
plot(smooth(mtss,5));
hold on;
plot(smooth(tss(69,:),5));
set(gca,'XTick',[1, 125,250,375,501]);
set(gca,'XTickLabel',{'-2.5kb','-1.25kb' , 'TSS', '+1.25kb' '2.5kb'})
ylabel('H3K4me3 density');%
title('WBC')
xlim([1,501]);
legend('pooled single cells','only one single cell');
