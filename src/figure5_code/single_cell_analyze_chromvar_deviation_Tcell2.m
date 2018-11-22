function [] = single_cell_analyze_chromvar_deviation_Tcell2(args1)

    a = readtable('../../data/input/Figure5/SRR161517_naive_Th1_Th2_mRNA_rep1_2.txt','Delimiter', '\t');
    tcell_gene = table2array(a(:,1));
    tcell_rpkm = table2array(a(:,2:7));
    import bioma.data.DataMatrix     
    DM0 = DataMatrix(log2(tcell_rpkm(:,1:2)+1));
    DM1 = DataMatrix(log2(tcell_rpkm(:,3:4)+1));
    DM2 = DataMatrix(log2(tcell_rpkm(:,5:6)+1));
    %DM0 = DataMatrix(tcell_rpkm(:,1:2));
    %DM1 = DataMatrix(tcell_rpkm(:,3:4));
    %DM2 = DataMatrix(tcell_rpkm(:,5:6));
    tt2 = [mean(tcell_rpkm(:,3:4)')',mean(tcell_rpkm(:,5:6)')',mean(tcell_rpkm(:,1:2)')'];
    %tt2 = [mean(tcell_rpkm(:,3:4)')',mean(tcell_rpkm(:,1:2)')'];

    [PValues1, TScores1] = mattest(DM1, DM0);
    [PValues2, TScores2] = mattest(DM2, DM0);
    rgene = tcell_gene;
    FValues1 = mafdr(PValues1);
    FValues2 = mafdr(PValues2);
    rgene = tcell_gene;
    %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    a = readtable('../../data/output/Figure5/chromvar/conet_het_scrna_sch3k4_rhigh_het_clow_het2_dev_tcell.txt');
    motif = table2array(a(:,1));
    dev_Th1_c1 = table2array(a(:,2:size(a,2)));
    a = readtable('../../data/output/Figure5/chromvar/conet_het_scrna_sch3k4_rhigh_het_chigh_het2_dev_tcell.txt');
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
    q1 = find(pval<0.05 & mean_diff'<-0.01); %0.045, 0.01
    q2 = find(pval<0.05 & mean_diff'>0.01); %0.045, 0.01
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
    [ia,ib,ic]=intersect(rgene,TF1);

    q = find(max(tt2(ib,:)')>1); % 7.5 for mono
    TF1 = ia(q);
    [ip,iq] = sort(tt2(ib(q),3)); %% change
    TF1 = TF1(iq);
    xbox1 = tt2(ib(q(iq)),:);


    [ia,ib,ic]=intersect(rgene,TF2);

    q = find(max(tt2(ib,:)')>1); % 7.5
    TF2 = ia(q);
    [ip,iq] = sort(tt2(ib(q),1)); %% change
    xbox2 = tt2(ib(q(iq)),:);
    TF2 = TF2(iq);

    xbox = [xbox1',xbox2']';


    if(args1==1)
        sTF1 = {'STAT6','BCL11B','PBX3'};
        sTF2 = {'PRDM1','E2F1'};

        [qs1,qs2, qs3] = intersect(TF1,sTF1);
        scatter(mean_diff,-log10(pval),'b','filled');
        hold on
        for i= 1:max(size(TF1(qs2)))
            ib = find(strcmp(TF_motif,TF1(qs2(i)))==1 & pval<0.05 & mean_diff'<-0.01);
            scatter(mean_diff(ib),-log10(pval(ib)),'r','filled');
            text(mean_diff(ib),-log10(pval(ib)),char(TF1(qs2(i))),'fontsize',12);
        end;
        [qs1,qs2, qs3] = intersect(TF2,sTF2);
        for i= 1:max(size(TF2(qs2)))
            ib2 = find(strcmp(TF_motif,TF2(qs2(i)))==1 & pval<0.05 & mean_diff'>0.01);
            scatter(mean_diff(ib2),-log10(pval(ib2)),'g','filled');
            text(mean_diff(ib2),-log10(pval(ib2)),char(TF2(qs2(i))),'fontsize',12);
        end    

        xlabel('Mean diff. in H3K4me3 levels associated with TFs');
        ylabel('-log10(Pvalue)');
        title('scChIC')
    end;

    if(args1==2)
        subplot(1,2,1);
        bar(xbox1(:,[1,3]))
        set(gca,'XTick',1:1:max(size(TF1)),'XTicklabel',TF1);
        set(gca,'XTickLabelRotation',45)
        ylabel('Gene expression');
        legend( 'Th1 cell','Naive T cell')

        subplot(1,2,2);
        bar(xbox2(:,[1,3]))
        set(gca,'XTick',1:1:max(size(TF2)),'XTicklabel',TF2);
        set(gca,'XTickLabelRotation',45)
        ylabel('Gene expression');
        legend( 'Th1 cell','Naive T cell')

        %imagesc(xbox);
        %set(gca,'XTick',[1,2],'XTicklabel',{'Th1','naive'},'YTick',1:1:max(size(TF1))+max(size(TF2)),'YTicklabel',[TF1',TF2'])
        %colorbar
    end;    
end
