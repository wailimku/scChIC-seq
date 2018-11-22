function [] = single_cell_analyze_chromvar_deviation_mono2(args1)


    a = readtable('../../data/input/Figure5/GSE36952_RPKM_3.txt','Delimiter', '\t');
    macro_gene = table2array(a(:,1));
    macro_rpkm = table2array(a(:,2:7));
    a = readtable('../../data/input/Figure5/GSM125682_combine3.txt','Delimiter', '\t');
    mono_gene = table2array(a(:,1));
    mono_rpkm = table2array(a(:,2:3));
    [ia,ib,ic] = intersect(mono_gene,macro_gene);
    mono_gene = ia;
    mono2_rpkm = mono_rpkm(ib,:);
    m1_rpkm = macro_rpkm(ic,1:3);
    m2_rpkm = macro_rpkm(ic,4:6);

    tt = [m1_rpkm, m2_rpkm, mono2_rpkm];
    tt2 = [mean(tt(:,1:3)')',mean(tt(:,4:6)')',mean(tt(:,7:8)')'];

    tt =log2(tt+1);


    import bioma.data.DataMatrix     
    DM0 = DataMatrix(tt(:,7:8));
    DM1 = DataMatrix(tt(:,1:3));
    DM2 = DataMatrix(tt(:,4:6));
    [PValues1, TScores1] = mattest(DM1, DM0);
    [PValues2, TScores2] = mattest(DM2, DM0);
    FValues1 = mafdr(PValues1);
    FValues2 = mafdr(PValues2);
    %=============================================================

    a = readtable('../../data/output/Figure5/chromvar/conet_het_scrna_sch3k4_rhigh_het_clow_het2_dev_mono.txt');
    motif = table2array(a(:,1));
    dev_Th1_c1 = table2array(a(:,2:size(a,2)));
    a = readtable('../../data/output/Figure5/chromvar/conet_het_scrna_sch3k4_rhigh_het_chigh_het2_dev_mono.txt');
    motif = table2array(a(:,1));
    dev_Th1_c2 = table2array(a(:,2:size(a,2)));




    TF_motif = cell([max(size(motif)),1]);
    for i = 1:max(size(motif))
        a = strsplit(char(motif(i)),'_');
        TF_motif(i) = cellstr(char(a(3)));
    end;    


    [qx,qy] = find(dev_Th1_c1==1000);
    %dev_Th1_c1(:,unique(qy))=[];
    dev_Th1_c1(unique(qx),:)=[];
    dev_Th1_c2(unique(qx),:)=[];
    motif(unique(qx))=[];
    TF_motif(unique(qx))=[];

    [qx,qy]= find(dev_Th1_c2==1000);
    %dev_Th1_c2(:,unique(qy))=[];
    %dev_Th1_c1(:,unique(qy))=[];

    dev_Th1_c2(unique(qx),:)=[];
    dev_Th1_c1(unique(qx),:)=[];
    motif(unique(qx))=[];
    TF_motif(unique(qx))=[];


    pval = zeros(max(size(motif)), 1);
    mean_diff = mean(dev_Th1_c1')-mean(dev_Th1_c2');
    for i = 1:max(size(motif))
        [h1,p1] = ttest2(dev_Th1_c1(i,:),dev_Th1_c2(i,:),'Vartype','equal','Tail','right');
        [h2,p2] = ttest2(dev_Th1_c1(i,:),dev_Th1_c2(i,:),'Vartype','equal','Tail','left');
        if(mean_diff(i)>0)
            pval(i) = p1;
        elseif(mean_diff(i)<0) 
            pval(i) = p2;
        end;    
    end;
    afdr = mafdr(pval);

    q1 = find(pval<0.05 & mean_diff'<-0.08); %0.05 -0.08
    q2 = find(pval<0.05 & mean_diff'>0.08);%0.05 0.08
    TF1 = unique(TF_motif(q1));
    TF2 = unique(TF_motif(q2));
    mm1 = motif(q1);
    mm2 = motif(q2);
    TF11 = setdiff(TF1,TF2); 
    TF22 = setdiff(TF2,TF1); 
    TF1 = TF11;
    TF2 = TF22;
    %TF1 = TF22;
    %TF2 = TF11;
    %[ia,ib,ic]=intersect(rgene,TF1);
    [ia,ib,ic]=intersect(mono_gene,TF1);

    q = find(max(tt2(ib,:)')>1); % 7.5 for mono
    TF1 = ia(q);
    [ip,iq] = sort(tt2(ib(q),1)); %% change
    TF1 = TF1(iq);
    %xbox1 = exp(log(2)*tt2(ib(q(iq)),:));
    xbox1 = tt2(ib(q(iq)),:);


    %[ia,ib,ic]=intersect(rgene,TF2);
    [ia,ib,ic]=intersect(mono_gene,TF2);

    q = find(max(tt2(ib,:)')>1); % 7.5
    TF2 = ia(q);
    [ip,iq] = sort(tt2(ib(q),3)); %% change
    %xbox2 = exp(log(2)*tt2(ib(q(iq)),:));
    xbox2 = tt2(ib(q(iq)),:);
    TF2 = TF2(iq);
    xbox = [xbox1',xbox2']';


    sTF1 = {'ENO1','ATF3','BACH1','MXI1'};
    sTF2 = {'FLI1','IKZF1','LYL1','STAT5A','CREB5','JDP2'};


    if(args1==1)
        [qs1,qs2, qs3] = intersect(TF1,sTF1);
        scatter(mean_diff,-log10(pval),'b','filled');
        hold on
        for i= 1:max(size(TF1(qs2)))
            ib = find(strcmp(TF_motif,TF1(qs2(i)))==1 & pval<0.05 & mean_diff'<-0.08);
            scatter(mean_diff(ib),-log10(pval(ib)),'r','filled');
            text(mean_diff(ib),-log10(pval(ib)),char(TF1(qs2(i))),'fontsize',12);
        end;
        [qs1,qs2, qs3] = intersect(TF2,sTF2);
        for i= 1:max(size(TF2(qs2)))
            ib2 = find(strcmp(TF_motif,TF2(qs2(i)))==1 & pval<0.05 & mean_diff'>0.08);
            scatter(mean_diff(ib2),-log10(pval(ib2)),'g','filled');
            text(mean_diff(ib2),-log10(pval(ib2)),char(TF2(qs2(i))),'fontsize',12);
        end    

        xlabel('Mean diff. in H3K4me3 levels associated with TFs');
        ylabel('-log10(Pvalue)');
        title('Monocyte');

    elseif(args1==2)
	subplot(1,2,1)
        bar(xbox1(:,[1,3]))
        set(gca,'XTick',1:1:max(size(TF1)),'XTicklabel',TF1);
        set(gca,'XTickLabelRotation',45)
        ylabel('Gene expression');
        legend( 'M1 Macrophage','Monocyte')
        %legend( 'Th1 cell','Naive T cell')

	subplot(1,2,2)
        bar(xbox2(:,[1,3]))
        set(gca,'XTick',1:1:max(size(TF2)),'XTicklabel',TF2);
        set(gca,'XTickLabelRotation',45)
        ylabel('Gene expression');
        legend( 'M1 Macrophage','Monocyte')
        %legend( 'Th1 cell','Naive T cell')
    end;

end
