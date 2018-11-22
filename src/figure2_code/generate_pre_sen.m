%a = readtable('../../data/kuw/biocore/wlku/Keji/scwbc_rc_52798_242_mat.txt','ReadVariableNames',0,'delimiter',',');
a = readtable('../../data/temp/Figure2/scwbc_rc_52798_242_mat.txt','ReadVariableNames',0,'delimiter',',');

b1 = table2array(a(:,2:243));
a2 = readtable('../../data/temp/Figure2/ff.txt','ReadVariableNames',0);
b2= table2array(a2);

ff = b2;
rc = b1;

c = a(:,1);
writetable(c,'../../data/temp/Figure2/test_peakname.txt','WriteVariableNames',0);
%========================================
%a = load('/data/kuw/biocore/wlku/Keji/scimpute/scwbc_52798_242.mat');
%ff = a.ff;
%rc = a.rc;
%========================================

pre2 = sum(rc)./ff';
sen2 = zeros(242,1);
for i = 1:242
	q = find(rc(:,i)>0);
	sen2(i) = max(size(q))*min(size(q))/52798;
end;

peak_len = 3000*ones(52798,1);
eff_genome = 3088286401*0.6;
nn = floor(eff_genome/mean(peak_len));
pp_peak = peak_len/eff_genome;
aa = zeros(52798,242);
for i = 1:242
    rr = ceil((rand(ff(i),1)*nn));
    q = find(rr<=52798);
    aa(rr(q),i) = aa(rr(q),i)+1;
end;

ran_pre = sum(aa)./ff';
ran_sen = zeros(242,1);
for i = 1:242
    q = find(aa(:,i)>0);
    ran_sen(i) = max(size(q))*min(size(q))/52798;
end;

[ia,ib] = sort(sen2,'descend');

pre40 = pre2(ib(1:24));
sen40 = sen2(ib(1:24));

gp1 = zeros(1,24+242+242);
gp1(1:24)=1;
gp1(25:242+24)=2;
gp1(243+24:24+242+242)=3;

%ffig3 = figure;

subplot(1,2,1);
boxplot([pre40, pre2, ran_pre],gp1);
ylabel('precision','fontsize',14);
set(gca,'XTick',[1,2,3],'XTicklabel',{'top10%', 'All','Random'});
subplot(1,2,2);
boxplot([sen40', sen2', ran_sen'],gp1);
ylabel('sensitivity','fontsize',14);
set(gca,'XTick',[1,2,3],'XTicklabel',{'top10%', 'All','Random'});


%saveas(ffig3,'../Figure/Figure2D_preandsen','pdf');
