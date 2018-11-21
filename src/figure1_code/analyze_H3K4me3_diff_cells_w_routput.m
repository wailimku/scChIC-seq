

path1 = '../../data/output/Figure1/';
path2 = '../../data/output/Figure1/';

a = readtable(strcat(char(path1),'bulk_diffcell_readcount_at_compeak.txt'));
a0 = table2array(a(:,2:4));

b = readtable(strcat(char(path1),'ChIC_diffcell_readcount_at_compeak.txt'));
b0 = table2array(b(:,2:4));


for i = 1:3
    a1(:,i) = a0(:,i)*1000000/sum(a0(:,i));
    b1(:,i) = b0(:,i)*1000000/sum(b0(:,i));
end;    
    

th2_rpkm = quantilenorm(log2(a1+1));
pth2_rpkm = quantilenorm(log2(b1+1));


ss = std(th2_rpkm(:,[1, 2, 3])');
qx = find(ss<=1.5);
th2_rpkm(qx,:)=[];
pth2_rpkm(qx,:)=[];
th2_zs = th2_rpkm;
pth2_zs = pth2_rpkm;

[ia,ib] = kmeans(th2_zs(:,[1, 2, 3]),6);


num_clus = 0;
for i =1:6
    q = find(ia==i);
    num_clus  = num_clus + floor(max(size(q))/200);
end;    
clus_th2_zs = zeros(num_clus,3);
pclus_th2_zs = zeros(num_clus,3);

layer1 = 100;

m=1;
for i = 1:6
    q = find(ia==i);
    qq = randperm(max(size(q)));
    n = 1;
    for j = 1:floor(max(size(q))/layer1)-1
        clus_th2_zs(m,:) = mean(th2_zs(q(qq(n*layer1+1:(n+1)*layer1)),[1,2,3]));
        pclus_th2_zs(m,:) = mean(pth2_zs(q(qq(n*layer1+1:(n+1)*layer1)),[1,2,3]));
        m = m+1;
        n = n+1;
    end;
     clus_th2_zs(m,:) = mean(th2_zs(q(qq(n*layer1+1:max(size(q)))),[1,2,3]));
     pclus_th2_zs(m,:) = mean(pth2_zs(q(qq(n*layer1+1:max(size(q)))),[1,2,3]));     
     m=m+1;
end;

subplot(1,2,1);
imagesc(pclus_th2_zs(:,[1,2,3]));
addpath('../../');
AdvancedColormap('bwr');
set(gca, 'XTick', 1:3);
set(gca, 'XTickLabel', {'3T3','ESC','naive'});
set(gca, 'YTickLabel', {});
colorbar;
title('3000 cells scChIC-seq');
caxis([prctile(pclus_th2_zs(:),5),prctile(pclus_th2_zs(:),95)]);

subplot(1,2,2);
imagesc(clus_th2_zs(:,[1,2,3]));
AdvancedColormap('bwr');
colorbar;
set(gca, 'XTick', 1:3);
set(gca, 'XTickLabel', {'3T3','ESC','naive'});
set(gca, 'YTickLabel', {});
title('Bulk cells ChIP-seq');
caxis([prctile(clus_th2_zs(:),5),prctile(clus_th2_zs(:),95)]);

