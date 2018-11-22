function single_cell_analyze_impute_readcount2(args1)

    order = [2,3,1,4]; 
    chr_length =[248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468];
    qqq1_chr1 = [1226, 1227, 1228, 1229, 1230, 1231, 1232, 1233, 1234, 1235, 1236, 1237, 1238, 1239, 1240, 1241, 1242, 1243, 1244, 1245, 1246, 1247, 1248, 1249, 1250, 1251, 1252, 1253, 1254, 1255, 1256, 1257, 1258, 1259, 1260, 1261, 1262, 1263, 1264, 1265, 1266, 1267, 1268, 1269, 1270, 1271, 1272, 1273, 1274, 1275, 1276, 1277, 1278, 1279, 1280, 1281, 1282, 1283, 1284, 1285, 1286, 1287, 1288, 1289, 1290, 1291, 1292, 1293, 1294, 1295, 1296, 1297, 1298, 1299, 1300, 1301, 1302, 1303, 1304, 1305, 1306, 1307, 1308, 1309, 1310, 1311, 1312, 1313, 1314, 1315, 1316, 1317, 1318, 1319, 1320, 1321, 1322, 1323, 1324, 1325, 1326, 1327, 1328, 1329, 1330, 1331, 1332, 1333, 1334, 1335, 1336, 1337, 1338, 1339, 1340, 1341, 1342, 1343, 1344, 1345, 1346, 1347, 1348, 1349, 1350, 1351, 1352, 1353, 1354, 1355, 1356, 1357, 1358, 1359, 1360, 1361, 1362, 1363, 1364, 1365, 1366, 1367, 1368, 1369, 1370, 1371, 1372, 1373, 1374, 1375, 1376, 1377, 1378, 1379, 1380, 1381, 1382, 1383, 1384, 1385, 1386, 1387, 1388, 1389, 1390, 1391, 1392, 1393, 1394, 1395, 1396, 1397, 1398, 1399, 1400, 1401, 1402, 1403, 1404, 1405, 1406, 1407, 1408, 1409, 1410, 1411, 1412, 1413, 1414, 1415];
    qqq1_chr7 =[592, 593, 594, 595, 596, 597, 598, 599, 600];
    qqq1_chr12 = [360,361,362,363,364,365,366,367,368];
    %bin = 2500:100000:248956422;
    bin = 2500:100000:chr_length(12);
    chr_rna = zeros(max(size(bin)),4);

    path = '../../data/input/Figure4/';

    %path = '/data/kuw/biocore/wlku/Keji/encode_wbc_rna/';
    file = {'GSM1256828_T_01_RPKM', 'GSM1256829_T_02_RPKM', 'GSM1256812_B_01_RPKM', 'GSM1256813_B_02_RPKM', 'GSM1256814_B_03_RPKM', 'GSM1256815_B_04_RPKM', 'GSM1256816_B_05_RPKM', 'GSM1256822_M_01_RPKM', 'GSM1256823_M_02_RPKM', 'GSM1256824_M_03_RPKM', 'GSM1256825_M_04_RPKM', 'GSM1256826_M_05_RPKM', 'GSM1982286_Sample_H1176_Rn372_HD_PBMC_NK_cell.RPKMforgenes', 'GSM1982287_Sample_H1176_Rn373_HD_PBMC_NK_cell.RPKMforgenes', 'GSM1982293_Sample_H1176_Rn379_HD_NK_cell.RPKMforgenes'};

    th_num = 21893;%
    b_num = 21893;
    mono_num = 21893;
    nk_num = 21850-4;

    th_rpkm = zeros(th_num, 2);
    th_gene = cell([th_num, 1]);
    th_chr = cell([th_num, 1]);
    th_pos = zeros(th_num, 1);

    b_rpkm = zeros(b_num, 5);
    b_gene = cell([b_num, 1]);
    b_chr = cell([b_num, 1]);
    b_pos = zeros(b_num, 1);


    mono_rpkm = zeros(mono_num, 5);
    mono_gene = cell([mono_num, 1]);
    mono_chr = cell([mono_num, 1]);
    mono_pos = zeros(mono_num, 1);

    nk_rpkm = zeros(nk_num, 3);
    nk_gene = cell([nk_num, 1]);
    nk_chr = cell([nk_num, 1]);
    nk_pos = zeros(nk_num, 1);
    %=========================================================================
    %=========================================================================
    for i = 1:2
        fp = fopen(strcat(char(path), char(file(i))), 'r');
        for j = 1:th_num
            a = fscanf(fp, '%s', 1);
            th_chr(j) = cellstr(a);    
            a = fscanf(fp, '%d', 1);
            b = fscanf(fp, '%d', 1);
            th_pos(j) = floor((a+b)/2);                
            a = fscanf(fp, '%s', 1);
            th_gene(j) = cellstr(a);
            th_rpkm(j,i) = fscanf(fp, '%f', 1);
            a = fscanf(fp, '%f', 1);      
        end;
        fclose(fp);
    end;    

    %for i = 1:2
    %    fp = fopen(strcat(char(path), char(file(i))), 'r');
    %    for j = 1:th_num
    %        a = fscanf(fp, '%s', 1);
    %        th_chr(j) = cellstr(a); 
    %        a = fscanf(fp, '%s', 1);
    %        th_gene(j) = cellstr(a);
    %        a = fscanf(fp, '%d', 1);
    %        b = fscanf(fp, '%d', 1);
    %        th_pos(j) = floor((a+b)/2);                
    %        th_rpkm(j,i) = fscanf(fp, '%f', 1);
    %        a = fscanf(fp, '%f', 1);      
    %    end;
    %    fclose(fp);
    %end;   


    for kkk = 3:7
        i = kkk-2;
        fp = fopen(strcat(char(path), char(file(kkk))), 'r');
        for j = 1:b_num
            a = fscanf(fp, '%s', 1);
            b_chr(j) = cellstr(a);    
            a = fscanf(fp, '%d', 1);
            b = fscanf(fp, '%d', 1);
            b_pos(j) = floor((a+b)/2);                
            a = fscanf(fp, '%s', 1);
            b_gene(j) = cellstr(a);
            b_rpkm(j,i) = fscanf(fp, '%f', 1);
            a = fscanf(fp, '%f', 1);        
        end;
        fclose(fp);    
    end;     



    for kkk = 8:12
        i = kkk-7;
        fp = fopen(strcat(char(path), char(file(kkk))), 'r');
        for j = 1:mono_num
            a = fscanf(fp, '%s', 1);
            mono_chr(j) = cellstr(a);    
            a = fscanf(fp, '%d', 1);
            b = fscanf(fp, '%d', 1);
            mono_pos(j) = floor((a+b)/2);                
            a = fscanf(fp, '%s', 1);
            mono_gene(j) = cellstr(a);
            mono_rpkm(j,i) = fscanf(fp, '%f', 1);
            a = fscanf(fp, '%f', 1);
        end;
        fclose(fp);    
    end;    

    for kkk = 13:15
        i = kkk-12;
        fp = fopen(strcat(char(path), char(file(kkk))), 'r');
        a = fgetl(fp);
        a = fgetl(fp);
        a = fgetl(fp);
        a = fgetl(fp);
        for j = 1:nk_num
            a = fscanf(fp, '%s', 1);
            nk_gene(j) = cellstr(a);  
            a = fscanf(fp, '%s', 1);
            nk_rpkm(j,i) = fscanf(fp, '%f', 1);        
            a = fscanf(fp, '%d', 1);        
        end;
        fclose(fp);
    end;    
    %=========================================================================
    %=========================================================================

    fb_rpkm = mean(b_rpkm,2);
    %fb_rpkm = b_rpkm(:,1);

    fmono_rpkm = mean(mono_rpkm,2);
    %fmono_rpkm = mono_rpkm(:,1);
    fth_rpkm = mean(th_rpkm,2);
    %fth_rpkm = th_rpkm(:,1);
    fnk_rpkm = mean(nk_rpkm(:,[1,2]),2);
    %fnk_rpkm = nk_rpkm(:,1);


    [ia,ib,ic] = intersect(mono_gene, nk_gene);
    fb_rpkm = fb_rpkm(ib);
    fmono_rpkm = fmono_rpkm(ib);
    fth_rpkm = fth_rpkm(ib);
    fnk_rpkm = fnk_rpkm(ic);

    nk_gene = nk_gene(ic);
    b_gene = b_gene(ib);
    th_gene = th_gene(ib);
    mono_gene = mono_gene(ib);


    %[ia,ib,ic] = intersect(mono_gene, th_gene);
    %fb_rpkm = fb_rpkm(ib);
    %fmono_rpkm = fmono_rpkm(ib);
    %fth_rpkm = fth_rpkm(ic);
    %fnk_rpkm = fnk_rpkm(ib);

    %nk_gene = nk_gene(ib);
    %b_gene = b_gene(ib);
    %th_gene = th_gene(ic);
    %mono_gene = mono_gene(ib);


    rpkmdata = [fb_rpkm, fmono_rpkm, fth_rpkm, fnk_rpkm];
    rpkmnormdata = rpkmdata;

    rpkmnormdata = quantilenorm(rpkmdata);


    a= readtable('../../data/input/Figure3/bulk_wbc_rc_52798_72_mat.txt');

    %a= readtable('/data/kuw/biocore/wlku/Keji/bulk_wbc_rc_52798_72_mat.txt');
    bulk_peakname= table2array(a(:,1));
    bulk_rc = table2array(a(:,2:73));

    bulk_cpm = bulk_rc;
    for i = 1:72
        bulk_cpm(:,i) = log2(bulk_rc(:,i)*1000000/sum(bulk_rc(:,i))+1);
    end;

    a= readtable('../../data/input/Figure3/bulk_wbc_file_name.txt');

    %a= readtable('/data/kuw/biocore/wlku/Keji/bulk_wbc_file_name.txt');
    bulkfile = table2array(a(:,2));

    wbc_cpm=bulk_cpm; 
    wbc_rc = bulk_rc;

    wbc_cpm_type = zeros(max(size(wbc_cpm)),12);
    %wbc_cpm_type(:,1) = mean(wbc_cpm(:,1:10),2);
    wbc_cpm_type(:,1) = mean(wbc_cpm(:,3:7),2);
    wbc_cpm_type(:,2) = mean(wbc_cpm(:,11:16),2);
    wbc_cpm_type(:,3) = mean(wbc_cpm(:,17:24),2);
    wbc_cpm_type(:,4) = mean(wbc_cpm(:,25:30),2);
    wbc_cpm_type(:,5) = mean(wbc_cpm(:,31:47),2);
    wbc_cpm_type(:,6) = wbc_cpm(:,48);
    wbc_cpm_type(:,7) = mean(wbc_cpm(:,49:50),2);
    %wbc_cpm_type(:,7) = mean(wbc_cpm(:,[17:24,35:30,31:47,49,50,62]),2);

    wbc_cpm_type(:,8) = mean(wbc_cpm(:,51:59),2);
    wbc_cpm_type(:,9) = mean(wbc_cpm(:,60:61),2);
    wbc_cpm_type(:,10) = wbc_cpm(:,62);
    wbc_cpm_type(:,11) = mean(wbc_cpm(:,63:64),2);
    wbc_cpm_type(:,12) = mean(wbc_cpm(:,65:72),2);

    wbc_rc_type = zeros(max(size(wbc_rc)),12);
    %wbc_rc_type(:,1) = mean(wbc_rc(:,1:10),2);
    wbc_rc_type(:,1) = mean(wbc_rc(:,3:7),2);
    wbc_rc_type(:,2) = mean(wbc_rc(:,11:16),2);
    wbc_rc_type(:,3) = mean(wbc_rc(:,17:24),2);
    wbc_rc_type(:,4) = mean(wbc_rc(:,25:30),2);
    wbc_rc_type(:,5) = mean(wbc_rc(:,31:47),2);
    wbc_rc_type(:,6) = wbc_rc(:,48);
    wbc_rc_type(:,7) = mean(wbc_rc(:,49:50),2);
    wbc_rc_type(:,8) = mean(wbc_rc(:,51:59),2);
    wbc_rc_type(:,9) = mean(wbc_rc(:,60:61),2);
    wbc_rc_type(:,10) = wbc_rc(:,62);
    wbc_rc_type(:,11) = mean(wbc_rc(:,63:64),2);
    wbc_rc_type(:,12) = mean(wbc_rc(:,65:72),2);

    wbc_rc_type(:,[3,4,6,7,8,9,12])=[];

    wbc_cpm_type(:,[3,4,6,7,8,9,12])=[];
    %================================%================================
    %========================================================================
    %========================================================================

    aaa1= readtable('../../data/input/Figure4/scimpute/scimpute_count.txt');
    aaa2= readtable('../../data/input/Figure4/peakname2.txt','Delimiter','\t','ReadVariableNames',0);
    sc3_cellclus = readtable('../../data/input/Figure4/cell_clus_sc3_12779_wimpute_nobatch_7a.txt');
    %aaa3= readtable('./data/input/Figure4/scwbc_rc_52798_242_mat.txt');

    %aaa1= readtable('../../data/input/Figure4/scimpute/scimpute_count.txt');
    %aaa2= readtable('../../data/input/Figure4/peakname2.txt','Delimiter','\t','ReadVariableNames',0);
    %sc3_cellclus = readtable('../../data/output/Figure3/SC3/cell_clus_sc3_wimpute_nobatch_7a.txt');


    %sclus = clus_idx(kkkk);
    rc = table2array(aaa1(:,2:243));
    peakname = table2array(aaa2);
    peakgene = peakname;
    dis = zeros(max(size(peakgene )),1);
    for i = 1:max(size(peakname))
        a = strsplit(char(peakname(i)),'_');
        peakgene(i) = cellstr(a(5));
        dis(i) = str2num(char(a(4)));
    end;

    %a = readtable('/data/kuw/biocore/wlku/Keji/scimpute/wbc_ident2.txt')
    %clus = table2array(a);
    rc(:,112)=[]; %12779 need to delete 112
    %rc(:,qdel)=[];
    clus = table2array(sc3_cellclus(:,2));

    cpm = rc;
    for i = 1:size(rc,2)
        cpm(:,i) = rc(:,i)*1000000/sum(rc(:,i));
    end;   

    q1 = find(clus==1); %% 1--T cell 7-- Mono
    cpm11 = cpm(:,q1);
    rc11 = rc(:,q1);
    q1 = find(clus==3); %% 1--T cell 7-- Mono
    cpm33 = cpm(:,q1);
    rc33 = rc(:,q1);
    q1 = find(clus==6); %% 1--T cell 7-- Mono
    cpm66 = cpm(:,q1);
    rc66 = rc(:,q1);
    q1 = find(clus==7); %% 1--T cell 7-- Mono
    cpm77 = cpm(:,q1);
    rc77 = rc(:,q1);

    q1 = find(sum(cpm11')==0 | sum(cpm33')==0 | sum(cpm66')==0 | sum(cpm77')==0);

    cpm11(q1,:)=[];
    cpm33(q1,:)=[];
    cpm66(q1,:)=[];
    cpm77(q1,:)=[];
    peakname2 = peakname;
    wbc_cpm_type2 = wbc_cpm_type;
    peakgene2 = peakgene;
    peakname2(q1)=[];
    wbc_cpm_type2(q1,:)=[];
    peakgene2(q1)=[];
    dis2 = dis;
    dis2(q1)=[];
    rc11(q1,:)=[];
    rc33(q1,:)=[];
    rc66(q1,:)=[];
    rc77(q1,:)=[];



    q = find(abs(dis2)>0);
    dis2(q)=[];
    cpm11(q,:)=[];
    cpm33(q,:)=[];
    cpm66(q,:)=[];
    cpm77(q,:)=[];

    peakgene2(q)=[];
    rc11(q,:)=[];
    rc33(q,:)=[];
    rc66(q,:)=[];
    rc77(q,:)=[];
    peakname2(q)=[];
    wbc_cpm_type2(q,:)=[];

    rr1 = corrcoef(log2(log2(cpm11'+1)+1));
    rr3 = corrcoef(log2(log2(cpm33'+1)+1));
    rr6 = corrcoef(log2(log2(cpm66'+1)+1));
    rr7 = corrcoef(log2(log2(cpm77'+1)+1));

    rrzs1 = rr1;
    rrzs3 = rr3;
    rrzs6 = rr6;
    rrzs7 = rr7;
    zs1 = zscore(rr1);
    zs3 = zscore(rr3);
    zs6 = zscore(rr6);
    zs7 = zscore(rr7);
    q1 = find(zs1<=0);
    q3 = find(zs3<=0);
    q6 = find(zs6<=0);
    q7 = find(zs7<=0);
    zs1(q1)=0;
    zs3(q3)=0;
    zs6(q6)=0;
    zs7(q7)=0;
	tic
    for i = 1:max(size(rr3))-1 %5732
        for j = i+1:max(size(rr3)) % 5733
            rrzs1(i,j) = sqrt(zs1(i,j).^2+ zs1(j,i).^2);
            rrzs3(i,j) = sqrt(zs3(i,j).^2+ zs3(j,i).^2);
            rrzs6(i,j) = sqrt(zs6(i,j).^2+ zs6(j,i).^2);
            rrzs7(i,j) = sqrt(zs7(i,j).^2+ zs7(j,i).^2);
            rrzs1(j,i) = rrzs1(i,j);
            rrzs3(j,i) = rrzs3(i,j);
            rrzs6(j,i) = rrzs6(i,j);
            rrzs7(j,i) = rrzs7(i,j);
        end;
	toc
    end;   
    rr1 = rrzs1;
    rr3 = rrzs3;
    rr6 = rrzs6;
    rr7 = rrzs7;


    rrc1 = rr1;
    q = find(rrc1<1.65);%% set 2.3
    rrc1(q)=0;
    q = find(rrc1>=1.65);%% set 2.3
    rrc1(q)=1; 

    rrc3 = rr3;
    q = find(rrc3<1.65);%% set 2.3
    rrc3(q)=0;
    q = find(rrc3>=1.65);%% set 2.3
    rrc3(q)=1;

    rrc6 = rr6;
    q = find(rrc6<1.65);%% set 2.3
    rrc6(q)=0;
    q = find(rrc6>=1.65);%% set 2.3
    rrc6(q)=1;

    rrc7 = rr7;
    q = find(rrc7<1.65);%% set 2.3
    rrc7(q)=0;
    q = find(rrc7>=1.65);%% set 2.3
    rrc7(q)=1;

    G1 = graph(rrc1);
    G3 = graph(rrc3);
    G6 = graph(rrc6);
    G7 = graph(rrc7);

    pg_rank1 = centrality(G1, 'Degree');
    pg_rank3 = centrality(G3, 'Degree');
    pg_rank6 = centrality(G6, 'Degree');
    pg_rank7 = centrality(G7, 'Degree');

    het11 = std(log2(cpm11'+1))./mean(log2(cpm11'+1));
    het33 = std(log2(cpm33'+1))./mean(log2(cpm33'+1));
    het66 = std(log2(cpm66'+1))./mean(log2(cpm66'+1));
    het77 = std(log2(cpm77'+1))./mean(log2(cpm77'+1));


    cometcut = 3;
    xxx = quantilenorm([pg_rank1,pg_rank3,pg_rank6,pg_rank7,het11',het33',het66',het77']);
    
    qk1 = find(xxx(:,1)./xxx(:,5)>cometcut);
    qk2 = find(xxx(:,5)./xxx(:,1)>cometcut);
    [ic2,id2] = sort(sum(cpm11'),'descend');
    [ie2,ig2] = sort(sum(cpm11'),'ascend');    
    qg1 = setdiff(qk1,id2(1:max(size(qk1))));
    qg2 = qk2;
    [ia,ib,ic] = intersect(th_gene,peakgene2(qg1));
    [ia2,ib2,ic2] = intersect(th_gene,peakgene2(qg2));
    gtu1 = peakgene2(qg1);
    gtd1 = peakgene2(qg2);

    qk1 = find(xxx(:,2)./xxx(:,6)>cometcut);
    qk2 = find(xxx(:,6)./xxx(:,2)>cometcut);
    [ic2,id2] = sort(sum(cpm33'),'descend');
    [ie2,ig2] = sort(sum(cpm33'),'ascend');    
    qg1 = setdiff(qk1,id2(1:max(size(qk1))));
    qg2 = qk2;
    [ia,ib,ic] = intersect(th_gene,peakgene2(qg1));
    [ia2,ib2,ic2] = intersect(th_gene,peakgene2(qg2));
    gtu3 = peakgene2(qg1);
    gtd3 = peakgene2(qg2);


    qk1 = find(xxx(:,3)./xxx(:,7)>cometcut);
    qk2 = find(xxx(:,7)./xxx(:,3)>cometcut);
    [ic2,id2] = sort(sum(cpm66'),'descend');
    [ie2,ig2] = sort(sum(cpm66'),'ascend');    
    qg1 = setdiff(qk1,id2(1:max(size(qk1))));
    qg2 = qk2;
    [ia,ib,ic] = intersect(th_gene,peakgene2(qg1));
    [ia2,ib2,ic2] = intersect(th_gene,peakgene2(qg2));
    gtu6 = peakgene2(qg1);
    gtd6 = peakgene2(qg2);

    qk1 = find(xxx(:,4)./xxx(:,8)>cometcut);
    qk2 = find(xxx(:,8)./xxx(:,4)>cometcut);
    [ic2,id2] = sort(sum(cpm77'),'descend');
    [ie2,ig2] = sort(sum(cpm77'),'ascend');    
    qg1 = setdiff(qk1,id2(1:max(size(qk1))));
    qg2 = qk2;
    [ia,ib,ic] = intersect(th_gene,peakgene2(qg1));
    [ia2,ib2,ic2] = intersect(th_gene,peakgene2(qg2));
    gtu7 = peakgene2(qg1);
    gtd7 = peakgene2(qg2);


    ppu1 = setdiff(setdiff(setdiff(setdiff(gtu1,gtu6),gtu7),gtu3),gtu7);
    ppd1 = setdiff(setdiff(setdiff(setdiff(gtd1,gtd6),gtu7),gtu3),gtd7);


    ppu7 = setdiff(setdiff(setdiff(gtu7,gtu6),gtu1),gtu3);
    ppd7 = setdiff(setdiff(setdiff(gtd7,gtd6),gtd1),gtu3);


    %ppu1 = setdiff(gtu1,gtu7);
    %ppd1 = setdiff(gtd1,gtd7);


    %ppu7 = setdiff(gtu7,gtu1);
    %ppd7 = setdiff(gtd7,gtd1);

    if(args1==1)
        [ia,ib,ic] = intersect(th_gene, ppu1); %% set ppu intersect
        [ia2,ib2,ic2] = intersect(th_gene,  ppd1);%%set ppd intersect
    elseif(args1==2)
        [ia,ib,ic] = intersect(th_gene, ppu7); %% set ppu intersect
        [ia2,ib2,ic2] = intersect(th_gene,  ppd7);%%set ppd intersect        
    end;    
        %[ia,ib,ic] = intersect(th_gene, setdiff(gtu1,gtu7)); %% set ppu intersect
    %[ia2,ib2,ic2] = intersect(th_gene,  setdiff(gtd1,gtd7));%%set ppd intersect


    cut1 = 6; %% set1
    cut2 = 6;
    Cpar = 0;
    q1 = find((rpkmnormdata(ib,1)>cut1));
    q2 = find((rpkmnormdata(ib2,1)>cut2));
    fprintf('%d\t%d\n',max(size(q1)),max(size(q2)));
    p1 = ranksum(log2(rpkmnormdata(ib(q1),1)+Cpar),log2(rpkmnormdata(ib2(q2),1)+Cpar));
    q1 = find((rpkmnormdata(ib,2)>cut1));
    q2 = find((rpkmnormdata(ib2,2)>cut2));
    fprintf('%d\t%d\n',max(size(q1)),max(size(q2)));
    p2 = ranksum(log2(rpkmnormdata(ib(q1),2)+Cpar),log2(rpkmnormdata(ib2(q2),2)+Cpar));
    q1 = find((rpkmnormdata(ib,3)>cut1));
    q2 = find((rpkmnormdata(ib2,3)>cut2));
    fprintf('%d\t%d\n',max(size(q1)),max(size(q2)));
    p3 = ranksum(log2(rpkmnormdata(ib(q1),3)+Cpar),log2(rpkmnormdata(ib2(q2),3)+Cpar));
    q1 = find((rpkmnormdata(ib,4)>cut1));
    q2 = find((rpkmnormdata(ib2,4)>cut2));
    fprintf('%d\t%d\n',max(size(q1)),max(size(q2)));
    p4 = ranksum(log2(rpkmnormdata(ib(q1),4)+Cpar),log2(rpkmnormdata(ib2(q2),4)+Cpar));
    pp=[p1,p2,p3,p4];


    celltype = {'B cell','Mono', 'T cell','NK cell'};
    for zzz = 1:4
        q1 = find((rpkmnormdata(ib,zzz)>cut1));
        q2 = find((rpkmnormdata(ib2,zzz)>cut1));  

        subplot(2,2,zzz);
        cdfplot(log2(rpkmnormdata(ib(q1),zzz)+Cpar));
        hold on;
        cdfplot(log2(rpkmnormdata(ib2(q2),zzz)+Cpar));
        title(strcat(char(celltype(zzz)),', pvalue =',num2str(round(pp(zzz),5,'significant'))));
        ylabel('cumulative distribution');
        
        
        xlabel('gene expression (log2 RPKM)');
        if(args1==1)
            legend('highly comet. peaks in Tcell','highly var. peaks in Tcell');
        elseif(args1==2)
            legend('highly comet. peaks in Mono','highly var. peaks in Mono');
        end;   
            
    end;


end
