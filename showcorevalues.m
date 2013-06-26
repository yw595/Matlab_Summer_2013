[excnumarray exctextarray raw]=xlsread('Supp Table 3 A community-driven global reconstruction of human metabolism 95.xls');
[height width]=size(excnumarray);
subexcnumarray=excnumarray(8:98,8:width);
smallestvalues=[];
allprctiles=[];
for i=2:2:120
    corevalues=subexcnumarray(:,i);
    sortedabscorevalues=sort(abs(corevalues));
    %disp(prctile(sortedabscorevalues,0:5:100));
    allprctiles(1:21,i/2)=prctile(sortedabscorevalues,0:5:100);
    %disp(sortedabscorevalues(6:10));
    %disp(median(sortedabscorevalues));
end
allprctiles
for i=1:size(allprctiles,1)
    median(allprctiles(i,:))
end