[excnumarray exctextarray raw]=xlsread('../Supp Table 3 A community-driven global reconstruction of human metabolism 95.xls');
[height width]=size(excnumarray);
subexcnumarray=excnumarray(8:98,8:width);
jainmetsarray=exctextarray(10:100,1);
for i=1:91
    disp([jainmetsarray{i} ' ' num2str(mean(subexcnumarray(i,:))*1/300)]);
end