[excnumarray exctextarray raw]=xlsread('../Supp Table 3 A community-driven global reconstruction of human metabolism 95.xls');
[height width]=size(excnumarray);
subexcnumarray=excnumarray(8:98,8:width);
celllinesarray=exctextarray(9,10:2:128);
outputFI=fopen('../comparezeros.txt','w');

metsarray={};
allinputvals1=[];
allinputvals2=[];
metstosigns=containers.Map;
for i=1:length(celllinesarray)
    if(~strcmp(celllinesarray{i},'MDA-MB-468') && ~strcmp(celllinesarray{i},'RXF 393'))
        expressionFile=strrep(celllinesarray{i},'(','_');
        expressionFile=strrep(expressionFile,')','_');
        expressionFile=strrep(expressionFile,' ','_');
        expressionFile=strrep(expressionFile,'/','_');
        expressionFile=strrep(expressionFile,'-','_');
        inputfile1=['../eMOMACorroutconstrainedless/' expressionFile 'out'];
        inputFI1=fopen(inputfile1,'r');
        inputfile2=['../eMOMACorroutconstrainedextra/' expressionFile 'out'];
        inputFI2=fopen(inputfile2,'r');
        
        line1=fgetl(inputFI1);
        line2=fgetl(inputFI2);
        inputvals1=[];
        inputvals2=[];
        flagexfluxes=0;
        while line1~=-1
            if(flagexfluxes)
                startindex=regexp(line1,'(\-|\d|\.)+$');
                if(length(metsarray)<91)
                    metsarray{end+1}=line1(1:startindex-2);
                end
                inputval1=str2num(line1(startindex:length(line1)));
                if(isKey(metstosigns,line1(1:startindex-2)))
                    valueArray=metstosigns(line1(1:startindex-2));
                    if(sign(inputval1)==1)
                        valueArray(1)=valueArray(1)+1;
                    elseif(sign(inputval1)==0)
                        valueArray(2)=valueArray(2)+1;
                    else
                        valueArray(3)=valueArray(3)+1;
                    end
                    metstosigns(line1(1:startindex-2))=valueArray;
                else
                    metstosigns(line1(1:startindex-2))=zeros(1,6);
                end
                inputvals1=[inputvals1; inputval1];            
            end
            if(~isempty(regexp(line1,'v_solex')))
                flagexfluxes=1;
            end
            line1=fgetl(inputFI1);
        end
        flagexfluxes=0;
        while line2~=-1
            if(flagexfluxes)
                startindex=regexp(line2,'(\-|\d|\.)+$');
                inputval2=str2num(line2(startindex:length(line2)));
                if(isKey(metstosigns,line2(1:startindex-2)))
                    valueArray=metstosigns(line2(1:startindex-2));
                    if(sign(inputval2)==1)
                        valueArray(4)=valueArray(4)+1;
                    elseif(sign(inputval2)==0)
                        valueArray(5)=valueArray(5)+1;
                    else
                        valueArray(6)=valueArray(6)+1;
                    end
                    metstosigns(line2(1:startindex-2))=valueArray;
                else
                    metstosigns(line2(1:startindex-2))=zeros(1,6);
                end
                inputvals2=[inputvals2; inputval2];
            end
            if(~isempty(regexp(line2,'v_solex')))
                flagexfluxes=1;
            end
            line2=fgetl(inputFI2);
        end
        allinputvals1=[allinputvals1 inputvals1];
        allinputvals2=[allinputvals2 inputvals2];
    end
end

metskeys=keys(metstosigns);
firstmetskeys={'alanine','arginine','asparagine','aspartate','biotin','folate','glucose','glutamine','glutamate','glutathione oxidized','glycine',...
'leucine','isoleucine','lysine','pantothenate','phenylalanine','proline','serine','tryptophan','tyrosine','valine'};
secondmetskeys={};
for i=1:length(metskeys)
    if(sum(strcmp(metskeys{i},firstmetskeys))==0)
        secondmetskeys{end+1}=metskeys{i};
    end
end
valAverage11=0;
valAverage12=0;
valAverage13=0;
for i=1:length(firstmetskeys)
    valueArray=metstosigns(firstmetskeys{i});
    valAverage11=valAverage11+valueArray(1)-valueArray(4);
    valAverage12=valAverage12+valueArray(2)-valueArray(5);
    valAverage13=valAverage13+valueArray(3)-valueArray(6);
    disp([firstmetskeys{i} ' ' num2str(valueArray(1)-valueArray(4)) ' ' num2str(valueArray(2)-valueArray(5)) ' ' num2str(valueArray(3)-valueArray(6))]);
end
disp('BREAK');
valAverage21=0;
valAverage22=0;
valAverage23=0;
for i=1:length(secondmetskeys)
    valueArray=metstosigns(secondmetskeys{i});
    valAverage21=valAverage21+valueArray(1)-valueArray(4);
    valAverage22=valAverage22+valueArray(2)-valueArray(5);
    valAverage23=valAverage23+valueArray(3)-valueArray(6);
    disp([secondmetskeys{i} ' ' num2str(valueArray(1)-valueArray(4)) ' ' num2str(valueArray(2)-valueArray(5)) ' ' num2str(valueArray(3)-valueArray(6))]);
end
valAverage11=valAverage11/length(firstmetskeys)
valAverage12=valAverage12/length(firstmetskeys)
valAverage13=valAverage13/length(firstmetskeys)
valAverage21=valAverage21/length(secondmetskeys)
valAverage22=valAverage22/length(secondmetskeys)
valAverage23=valAverage23/length(secondmetskeys)
[allinputvalsheight allinputvalswidth]=size(allinputvals1);
for i=1:allinputvalsheight
    allinputvals1row=allinputvals1(i,:);
    allinputvals2row=allinputvals2(i,:);
    fprintf(outputFI,'%s\t',metsarray{i});
    for j=1:allinputvalswidth
        fprintf(outputFI,'%d %f %f\t', j,allinputvals1row(j),allinputvals2row(j));
    end
    fprintf(outputFI,'\n');
end