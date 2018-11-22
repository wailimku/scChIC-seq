function []= single_cell_analyze_corr_scRNA_in_concen_mat2(args1, args2, args3,  args4)

    a= readtable('../../data/input/Figure4/scimpute_Thscrna/scimpute_count.txt');
    %a= readtable('/data/kuw/biocore/wlku/Keji/scimpute_Thscrna/scimpute_count.txt');
    %a= readtable('/data/kuw/biocore/wlku/Keji/scimpute/Th_cell_seurat_norm_data.txt');
    Th_gene = table2array(a(:,1));
    Th_cnt = table2array(a(:,2:191));
    Th_exp = Th_cnt;
    for i = 1:190
        Th_exp(:,i) = log2(Th_exp(:,i)*1000000/sum(Th_exp(:,i))+1);
    end;
    Th_het = std(Th_exp')./mean(Th_exp');

    a= readtable('../../data/input/Figure4/scimpute_Bscrna/scimpute_count.txt');
    %a= readtable('/data/kuw/biocore/wlku/Keji/scimpute_Bscrna/scimpute_count.txt');
    %a= readtable('/data/kuw/biocore/wlku/Keji/scimpute/B_cell_seurat_norm_data.txt');
    B_gene = table2array(a(:,1));
    B_cnt = table2array(a(:,2:331));
    B_exp = B_cnt;

    for i = 1:330
        B_exp(:,i) = log2(B_exp(:,i)*1000000/sum(B_exp(:,i))+1);
    end;
    B_het = std(B_exp')./mean(B_exp');

    a= readtable('../../data/input/Figure4/scimpute_NKscrna/scimpute_count.txt');
    %a= readtable('/data/kuw/biocore/wlku/Keji/scimpute_NKscrna/scimpute_count.txt');
    %a= readtable('/data/kuw/biocore/wlku/Keji/scimpute/NK_cell_seurat_norm_data.txt');
    NK_gene = table2array(a(:,1));
    NK_cnt = table2array(a(:,2:769));
    NK_exp = NK_cnt;
    for i = 1:768
        NK_exp(:,i) = log2(NK_exp(:,i)*1000000/sum(NK_exp(:,i))+1);
    end;
    NK_het = std(NK_exp')./mean(NK_exp');

    a= readtable('../../data/input/Figure4/scimpute_Monoscrna/scimpute_count.txt');
    %a= readtable('/data/kuw/biocore/wlku/Keji/scimpute_Monoscrna/scimpute_count.txt');
    %a= readtable('/data/kuw/biocore/wlku/Keji/scimpute/Mono_seurat_norm_data.txt');
    Mono_gene = table2array(a(:,1));
    Mono_cnt = table2array(a(:,2:130));
    Mono_exp = Mono_cnt;
    for i = 1:129
        Mono_exp(:,i) = log2(Mono_exp(:,i)*1000000/sum(Mono_exp(:,i))+1);
    end;
    Mono_het = std(Mono_exp')./mean(Mono_exp');
    %================================%================================
    %========================================================================
    %========================================================================
    %========================================================================
    a= readtable('../../data/input/Figure4/scimpute/scimpute_count.txt');
    %a= readtable('/data/kuw/biocore/wlku/Keji/scimpute/scimpute_count.txt');
    %a= readtable('/data/kuw/biocore/wlku/Keji/scwbc_rc_52798_242_mat.txt');
    rc = table2array(a(:,2:243));
    a= readtable('../../data/input/Figure4/scimpute/peakname2.txt','Delimiter','\t','ReadVariableNames',0);
    %a= readtable('/data/kuw/biocore/wlku/Keji/scimpute/peakname2.txt','Delimiter','\t','ReadVariableNames',0);
    peakname = table2array(a);
    peakgene = peakname;
    dis = zeros(max(size(peakgene )),1);
    for i = 1:max(size(peakname))
        a = strsplit(char(peakname(i)),'_');
        peakgene(i) = cellstr(a(5));
        dis(i) = str2num(char(a(4)));
    end;

    %a = readtable('/data/kuw/biocore/wlku/Keji/scimpute/wbc_ident2.txt')
    %clus = table2array(a);

    rc(:,112)=[]; % 12779**********
    %rc(:,qdel)=[];
    %sc3_cellclus = readtable('../../data/output/Figure3/SC3/cell_clus_sc3_wimpute_nobatch_8a.txt');
    sc3_cellclus = readtable('../../data/input/Figure4/cell_clus_sc3_12779_wimpute_nobatch_7a.txt');
    %sc3_cellclus = readtable('/data/kuw/biocore/wlku/Keji/wbc_285_read1_mapped/SC3/cell_clus_sc3_12779_wimpute_nobatch_7a.txt');
    clus = table2array(sc3_cellclus(:,2));

    cpm = rc;
    for i = 1:size(rc,2)
        cpm(:,i) = rc(:,i)*1000000/sum(rc(:,i));
    end;    
    q = find(clus==args1);  %%  T cell ==1 Mono==7
    cpm33 = cpm(:,q);
    rc33 = rc(:,q);
    q1 = find(sum(cpm33')==0);
    cpm33(q1,:)=[];
    peakname2 = peakname;
    peakgene2 = peakgene;
    peakname2(q1)=[];
    peakgene2(q1)=[];
    dis2 = dis;
    dis2(q1)=[];
    rc33(q1,:)=[];

    q = find(abs(dis2)>0);
    dis2(q)=[];
    cpm33(q,:)=[];
    peakgene2(q)=[];
    rc33(q,:)=[];
    peakname2(q)=[];

    %+==============================================================
    if(args1==1)
        [ia2,ib2,ic2] = intersect(Th_gene, peakgene2,'stable');
        scchip_mat = log2(cpm33(ic2,:)+1);
        scrna_mat = Th_exp(ib2,:);
    elseif(args1==7)
        [ia2,ib2,ic2] = intersect(Mono_gene, peakgene2,'stable');
        scchip_mat = log2(cpm33(ic2,:)+1);
        scrna_mat = Mono_exp(ib2,:);             
    end;    
    
    peakname3 = peakname2(ic2);
    q = find(sum(scrna_mat')==0|sum(scchip_mat')==0);
    scrna_mat(q,:)=[];
    scchip_mat(q,:)=[];
    peakname3(q)=[];
    peakgene3 = ia2;
    peakgene3(q)=[];


    rhet = std(scrna_mat')./mean(scrna_mat');
    chet = std(scchip_mat')./mean(scchip_mat');

    if(args2==1)
        q = find(chet<0);
        rhet(q)=[];
        chet(q)=[];
        xx = chet;
        q1 = find(xx<=prctile(xx,100/3));
        q2 = find(xx>prctile(xx,100/3) & xx<prctile(xx,200/3));
        q3 = find(xx>=prctile(xx,200/3));

        minsize = min(min(max(size(q1)),max(size(q2))),max(size(q3)));
        if(max(size(q1))>minsize)
           q1(minsize+1:max(size(q1)))=[];
        end;
        if(max(size(q2))>minsize)
           q2(minsize+1:max(size(q2)))=[];
        end;
        if(max(size(q3))>minsize)
           q3(minsize+1:max(size(q3)))=[];
        end;
        xx1 = [mean(chet(q1)),mean(chet(q2)),mean(chet(q3))];
        xx2 = [mean(rhet(q1)),mean(rhet(q2)),mean(rhet(q3))];


        %q1(2268) =[]; % mono
        %q3(2268) =[];% mono
        %q1(2319) =[]; % mono
        %q3(2319:2320) =[];% mono


        Y = zeros(3, max(size(q1)));
        Y(1,:) = rhet(q1);
        Y(2,:) = rhet(q2);
        Y(3,:) = rhet(q3);
        Y=Y';

        ranksum(rhet(q1), rhet(q2))
        ranksum(rhet(q3), rhet(q2))
        addpath('../../');
        violinplot(Y,{'bottom 1/3','middle 1/3','top 1/3'},'MedianColor',[0,0,0])
        ylim([0,17]);
        ylabel('Cofficient of variation(scRNA)');
        title('scChIC');
	ylim([args3,args4]);

    elseif(args2==2)


        rr1 = corrcoef(log2(scrna_mat'+1));
        rr3 = corrcoef(log2(scchip_mat'+1));

        %rr2 = corrcoef(scrna_mat');
        %rr4 = corrcoef(scchip_mat');

        rrzs1 = rr1;
        rrzs3 = rr3;
        zs1 = zscore(rr1);
        zs3 = zscore(rr3);
        q1 = find(zs1<=0);
        q3 = find(zs3<=0);
        zs1(q1)=0;
        zs3(q3)=0;
        tic
        for i = 1:max(size(rr3))-1 %5732
            for j = i+1:max(size(rr3)) % 5733
                rrzs1(i,j) = sqrt(zs1(i,j).^2+ zs1(j,i).^2);
                rrzs3(i,j) = sqrt(zs3(i,j).^2+ zs3(j,i).^2);
                rrzs1(j,i) = rrzs1(i,j);
                rrzs3(j,i) = rrzs3(i,j);
            end;
            toc;
        end;   

        rr1 = rrzs1; % scrna
        rr3 = rrzs3; % scChIC

        q = find(rr3<1.65);
        rr3(q)=[];
        rr1(q)=[];

        xx = rr3;
        q1 = find(xx<=prctile(xx,100/3));
        q2 = find(xx>prctile(xx,100/3) & xx<=prctile(xx,200/3));
        q3 = find(xx>=prctile(xx,200/3));

        minsize = min(min(max(size(q1)),max(size(q2))),max(size(q3)));

        if(max(size(q1))>minsize)
           q1(minsize+1:max(size(q1)))=[];
        end;
        if(max(size(q2))>minsize)
           q2(minsize+1:max(size(q2)))=[];
        end;
        if(max(size(q3))>minsize)
           q3(minsize+1:max(size(q3)))=[];
        end;

        xx1 = [mean(rr3(q1)),mean(rr3(q2)),mean(rr3(q3))];
        xx2 = [mean(rr1(q1)),mean(rr1(q2)),mean(rr1(q3))];




        Y = zeros(3, max(size(q1)));
        Y(1,:) = rr1(q1);
        Y(2,:) = rr1(q2);
        Y(3,:) = rr1(q3);
        Y=Y';
	addpath('../../')
        violinplot(Y,{'bottom 1/3','middle 1/3','top 1/3'},'MedianColor',[0,0,0])
        %ylim([0,17]);
        ylabel('Co-expression(scRNA)');
        title('scChIC (Co-methylation)')
	ylim([args3,args4]);

        ranksum(rr1(q1),rr1(q2))
        ranksum(rr1(q2),rr1(q3))
    end;   
end
