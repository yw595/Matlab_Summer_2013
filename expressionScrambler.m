function expressionScrambler(expressionFile, outputFile)
inputFI=fopen(expressionFile,'r');
%expressionFile
outputFI=fopen(outputFile,'w');
line=fgetl(inputFI);
line=fgetl(inputFI); %move to first line after headers
geneIDs={};
geneMeans={};
geneVars={};

while line~=-1
    [startIndex1,endIndex1]=regexp(line,'^(\w)+\t');
    geneID=strtrim(line(startIndex1:endIndex1));
    [startIndex2,endIndex2]=regexp(line,'\t(nan|((\d)*(\.)?(\d)+))\t');
    geneMean=strtrim(line(startIndex2:endIndex2));
    [startIndex3,endIndex3]=regexp(line,'\t(nan|((\d)*(\.)?(\d)+))$');
    geneVar=strtrim(line(startIndex3:endIndex3));
    geneIDs{end+1}=geneID;geneMeans{end+1}=geneMean;geneVars{end+1}=geneVar;
    %disp([geneID ' ' geneMean ' ' geneVar]);
    line=fgetl(inputFI);
end

randomInds=randperm(length(geneMeans));
scrambledgeneMeans=geneMeans(randomInds);
scrambledgeneVars=geneVars(randomInds);
fprintf(outputFI,'gene\tmean\tvar\n');
for i=1:length(geneIDs);
    fprintf(outputFI,[geneIDs{i} '\t' scrambledgeneMeans{i} '\t' scrambledgeneVars{i} '\n']);
end
fclose(outputFI);
end

