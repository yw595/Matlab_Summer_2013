[excnumarray exctextarray raw]=xlsread('../Supp Table 3 A community-driven global reconstruction of human metabolism 95.xls');
celllinesarray=exctextarray(9,10:2:128);

for i=1:length(celllinesarray)
    if(~strcmp(celllinesarray{i},'MDA-MB-468') && ~strcmp(celllinesarray{i},'RXF 393'))
        expressionFile=strrep(celllinesarray{i},'(','_');
        expressionFile=strrep(expressionFile,')','_');
        expressionFile=strrep(expressionFile,' ','_');
        expressionFile=strrep(expressionFile,'/','_');
        expressionFile=strrep(expressionFile,'-','_');
        coreFile = ['../NCI60exp/' expressionFile '.csv'];
        outFile = ['../NCI60expscrambled/' expressionFile '.csv'];
        expressionScrambler(coreFile,outFile);
    end
end