[excarray textarray raw]=xlsread('Supp Table 3 A community-driven global reconstruction of human metabolism 95.xls');
[height width]=size(excarray);
subexcarray=excarray(8:98,8:width);
fvaarray=excarray(8:98,[1 3]);
[spearmanrho spearmanpval]=corr(fvaarray,subexcarray,'type','Spearman');
[pearsonrho pearsonpval]=corr(fvaarray,subexcarray,'type','Pearson');
[kendallrho kendallpval]=corr(fvaarray,subexcarray,'type','Kendall');