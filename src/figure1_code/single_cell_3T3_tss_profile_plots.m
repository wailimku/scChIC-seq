kkk=2;

if(kkk==1)
    addpath('./data/output/Figure1/');
    a = readtable('3T3_lowcell_output.txt','ReadVariableNames',0);
    tss = zeros(5,501);
    for i = 1:5
        b = str2num(str2mat(cellstr(table2array(a(2:502,2+(i-1)*3)))));
        tss(i,:) = b';
    end;
elseif(kkk==2)
    a = readtable('./data/temp/Figure1/3T3_low_cell_output_scChIC.txt','ReadVariableNames',0);
    b = readtable('./data/input/Figure1/3T3_bulk_output.txt','ReadVariableNames',0);
    tss = zeros(5,501);
    c = str2num(str2mat(cellstr(table2array(b(2:502,2)))));
    tss(1,:) = c';
   for i = 1:4
        c = str2num(str2mat(cellstr(table2array(a(2:502,2+(i-1)*3)))));
        tss(i+1,:) = c';
    end;
end;

fig = figure;

plot(smooth(tss(2,:),10)/1000,'r');
hold on
plot(smooth(tss(3,:),10)/1000,'b');
plot(smooth(tss(4,:),10)/1000,'m');
plot(smooth(tss(5,:),10)/1000,'g');
plot(smooth(tss(1,:),10)/1000,'k');


xlim([1,501]);
legend('3000 cells','1000 cells','300 cells','100 cells','Bulk')

set(gca,'XTick',[1, 125,250,375, 501]);
set(gca,'XTickLabel',{'-5kb','-2.5kb' , 'TSS', '+2.5kb' '5kb'});
ylabel('Density');
xlim([1,501]);

%print(fig,'./Figures/Figure1/3T3_TSS_profiles.pdf','-dpdf')
