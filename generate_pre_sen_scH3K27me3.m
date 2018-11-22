%========================================
%a = load('/data/kuw/biocore/wlku/pipeline/run1070_KZ1655_KZ1656/temp_selbed/sel_preci_012_bed/scwbc_52798_242.mat');
%ff = a.ff;
%rc = a.rc;
%a = readtable('/data/kuw/biocore/wlku/pipeline/run1070_KZ1655_KZ1656/temp_selbed/sel_preci_012_bed/bulk_wbc_rc_36224_106_mat_2.txt','ReadVariableNames',0);
a = readtable('../../data/temp/Figure6/bulk_wbc_rc_36224_106_mat_2.txt','ReadVariableNames',0);
rc = table2array(a(:,2:107));
%a2 = readtable('/data/kuw/biocore/wlku/pipeline/run1070_KZ1655_KZ1656/temp_selbed/sel_preci_012_bed/sc_wbc_file_name_2.txt','ReadVariableNames',0,'delimiter',' ');
a2 = readtable('../../data/temp/Figure6/sc_wbc_file_name_2.txt','ReadVariableNames',0,'delimiter',' ');

b2 = table2array(a2(:,2));
%a3 = readtable('/data/kuw/biocore/wlku/pipeline/run1070_KZ1655_KZ1656/temp_selbed/sel_preci_012_bed/sel_filterlist_1000_08_3.txt','ReadVariableNames',0,'delimiter','\t');
a3 = readtable('../../data/temp/Figure6/sel_filterlist_106.txt','ReadVariableNames',0,'delimiter','\t');

ff = table2array(a3(:,1));
b3 = table2array(a3(:,2));

[ia,ib,ic] = intersect(b2,b3,'stable');
ff = ff(ic);
%========================================

pre2 = sum(rc)./ff';
sen2 = zeros(106,1);
for i = 1:106
	q = find(rc(:,i)>0);
	sen2(i) = max(size(q))*min(size(q))/36224;
end;

q2 = find(ff<1000);
pre3 = pre2;
sen3 = sen2;
ff2 = ff;

pre3(q2)=[];
sen3(q2)=[];
ff2(q2)=[];


peak_len = 4000*ones(36224,1);
eff_genome = 3088286401*0.8;
nn = floor(eff_genome/mean(peak_len));
pp_peak = peak_len/eff_genome;
aa = zeros(36224,84);
for i = 1:84
    rr = ceil((rand(ff2(i),1)*nn));
    q = find(rr<=36224);
    aa(rr(q),i) = aa(rr(q),i)+1;
end;

ran_pre = sum(aa)./ff2';
ran_sen = zeros(84,1);
for i = 1:84
    q = find(aa(:,i)>0);
    ran_sen(i) = max(size(q))*min(size(q))/36224;
end;


gp1 = zeros(1,168);
gp1(1:84)=1;
gp1(84:168)=2;


subplot(1,2,1);
boxplot([pre3, ran_pre],gp1);
ylabel('precision','fontsize',14);
set(gca,'XTick',[1,2],'XTicklabel',{'All','Random'});

subplot(1,2,2);
boxplot([sen3', ran_sen'],gp1);
ylabel('sensitivity','fontsize',14);
ylim([0,0.25])
set(gca,'XTick',[1,2],'XTicklabel',{ 'All','Random'});


%saveas(ffig3,'../Figure/Figure2D_preandsen2','pdf');