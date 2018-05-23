%========================================
a = load('../data/output/single_cell_rc_data.mat');
ff = a.ff;
sc_peak_id = a.sc_peak_id;
chr = a.chr;
peakss = a.peakss;
peakes = a.peakes;
rc = a.rc;
cpm = a.cpm;
ib_spec  =a.ib_spec;
ic_spec = a.ic_spec;
clear a;
%========================================
a = load('../data/input/H3K4me3_wbc_ChIP-seq_ENCODE2.mat');
wbc_id = a.wbc_id;
wbc_rc = a.wbc_rc;
wbc_cpm_type = a.wbc_cpm_type;
wbc_rc_type = a.wbc_rc_type;
wbc_type_name = a.wbc_type_name;
clear a;

pre2 = sum(rc)./ff';
sen2 = zeros(242,1);
for i = 1:242
	q = find(rc(:,i)>0);
	sen2(i) = max(size(q))*min(size(q))/41426;
end;

peak_len = peakes - peakss;
eff_genome = 3088286401*0.6;
nn = floor(eff_genome/mean(peak_len));
pp_peak = peak_len/eff_genome;
aa = zeros(41426,242);
for i = 1:242
    rr = ceil((rand(ff(i),1)*nn));
    q = find(rr<=41426);
    aa(rr(q),i) = aa(rr(q),i)+1;
end;

ran_pre = sum(aa)./ff';
ran_sen = zeros(242,1);
for i = 1:242
    q = find(aa(:,i)>0);
    ran_sen(i) = max(size(q))/41426;
end;

[ia,ib] = sort(sen2,'descend');

pre40 = pre2(ib(1:40));
sen40 = sen2(ib(1:40));

gp1 = zeros(1,242+242+40);
gp1(1:40)=1;
gp1(41:40+242)=2;
gp1(283:40+242+242)=3;

ffig3 = figure;

subplot(1,2,1);
boxplot([pre40, pre2, ran_pre],gp1);
ylabel('precision','fontsize',14);
subplot(1,2,2);
boxplot([sen40', sen2', ran_sen'],gp1);
ylabel('sensitivity','fontsize',14);

saveas(ffig3,'../Figure/Figure2D_preandsen','pdf');