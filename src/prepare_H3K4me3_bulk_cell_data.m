
%===========================================================
%
%===========================================================

path1 = '../data/input/';
a=load(strcat(char(path1),'H3K4me3_wbc_ChIP-seq_ENCODE.mat'));
wbc_id = a.wbc_id;
wbc_cpm_type = a.wbc_cpm_type;
wbc_rc_type = a.wbc_rc_type;
wbc_type_name = a.wbc_type_name;
clear a;
