function expressionModifier(expressionFile, geneID, factor)
inputFI=fopen(expressionFile,'r');
outputFile=strcat(expressionFile,'2');
outputFI=fopen(outputFile,'w');
line=fgetl(inputFI);
while line~=-1
    if(~isempty(regexp(line,['^',num2str(geneID),'\t'])))
        %a='HERE'
        [startIndex1,endIndex1]=regexp(line,'\t(.)*\t');
        mean=strtrim(line(startIndex1:endIndex1));
        %mean
        %size(mean)
        newmean=factor*str2num(line(startIndex1:endIndex1));
        [startIndex2,endIndex2]=regexp(line,'\t(\d)*.(\d)*$');
        sd=strtrim(line(startIndex2:endIndex2));
        %sd
        %size(sd)
        line=[num2str(geneID),'\t',num2str(newmean),'\t',sd];
    end
    fprintf(outputFI,strcat(line,'\n'));
    line=fgetl(inputFI);
end
end

