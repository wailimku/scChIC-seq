
addpath('./data/input/Figure1/');
a = readtable('3T3_ESC_T_homer_output.txt','ReadVariableNames',0);

tss = zeros(6,501);
for i = 1:6
    b = str2num(str2mat(cellstr(table2array(a(2:502,2+(i-1)*3)))));
    tss(i,:) = b';
end;

subplot(1,3,1)
plot(smooth(tss(4,:),10)/1000);
hold on
plot(smooth(tss(1,:),10)/1000);
xlim([1,501]);
legend('H3K4me3 Ab-MNase (3000 cells)','IgG Ab-Mnase (3000 cells)')
set(gca,'XTick',[1, 125,250,375, 501]);
set(gca,'XTickLabel',{'-5kb','-2.5kb' , 'TSS', '+2.5kb' '5kb'});
ylabel('Density');
title('3T3');
ylim([0,0.07])

subplot(1,3,2)
plot(smooth(tss(5,:),10)/1000);
hold on
plot(smooth(tss(2,:),10)/1000);
xlim([1,501]);
legend('H3K4me3 Ab-MNase (3000 cells)','IgG Ab-Mnase (3000 cells)')
set(gca,'XTick',[1, 125,250,375, 501]);
set(gca,'XTickLabel',{'-5kb','-2.5kb' , 'TSS', '+2.5kb' '5kb'});
ylabel('Density');
title('mESC');
ylim([0,0.07])

subplot(1,3,3)
plot(smooth(tss(6,:),10)/1000);
hold on
plot(smooth(tss(3,:),10)/1000);
xlim([1,501]);
legend('H3K4me3 Ab-MNase (3000 cells)','IgG Ab-Mnase (3000 cells)')
set(gca,'XTick',[1, 125,250,375, 501]);
set(gca,'XTickLabel',{'-5kb','-2.5kb' , 'TSS', '+2.5kb' '5kb'});
ylabel('Density');
title('naive T');
ylim([0,0.07])



xlim([1,501]);
legend('H3K4me3 Ab-MNase','IgG Ab-Mnase')

set(gca,'XTick',[1, 125,250,375, 501]);
%set(gca,'XTick',[1,100]);
set(gca,'XTickLabel',{'-5kb','-2.5kb' , 'TSS', '+2.5kb' '5kb'});
ylabel('Density');
xlim([1,501]);
saveas(gcf,'./Figures/Figure1/3T3_TSS_profiles.pdf');