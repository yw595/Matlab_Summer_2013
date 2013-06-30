[excnumarray exctextarray raw]=xlsread('../Supp Table 3 A community-driven global reconstruction of human metabolism 95.xls');
[height width]=size(excnumarray);
subexcnumarray=excnumarray(8:98,8:width);
smallestvalues=[];
allprctiles=[];
celllinesarray=exctextarray(9,10:2:128);
celllinesarray2={};
allnumzerodiffs=[];
allnumcorezero=[];
for i=1:length(celllinesarray)
    if(~strcmp(celllinesarray{i},'MDA-MB-468') && ~strcmp(celllinesarray{i},'RXF 393'))
        expressionFile=strrep(celllinesarray{i},'(','_');
        expressionFile=strrep(expressionFile,')','_');
        expressionFile=strrep(expressionFile,' ','_');
        expressionFile=strrep(expressionFile,'/','_');
        expressionFile=strrep(expressionFile,'-','_');
        celllinesarray2{end+1}=celllinesarray{i};
        
        inputfile=strcat(strcat('../eMOMACorrout2/',expressionFile),'out');
        inputFI=fopen(inputfile,'r');
        corevalues=subexcnumarray(:,i*2);
        abscorevalues=abs(corevalues);
        [sortedabscorevalues sortedabscoreinds]=sort(abs(corevalues));
        
        line=fgetl(inputFI);
        calculatedfluxes=[];
        flagexfluxes=0;
        linenum=0;
        numcalculatedzero=0;
        numcalculatedandcorezero=0;
        numcorezero=0;
        %j=1;
        while line~=-1
            if(flagexfluxes)
                startindex=regexp(line,'(\-|\d|\.)+$');
                linenum=linenum+1;
                if(abscorevalues(linenum)<=.0001)
                    numcorezero=numcorezero+1;
                end
                if(str2num(line(startindex:length(line)))==0)
                    numcalculatedzero=numcalculatedzero+1;
                    if(abscorevalues(linenum)<=.0001)
                        numcalculatedandcorezero=numcalculatedandcorezero+1;
                    end
                end
                calculatedfluxes=[calculatedfluxes; str2num(line(startindex:length(line)))];
            end
            if(~isempty(regexp(line,'v_solex')))
                flagexfluxes=1;
            end
            line=fgetl(inputFI);
        end
        %numcalculatedzero
        %numcorezero
        numcalculatedandcorezero;
        numcorezero*numcalculatedzero/91;
        allnumzerodiffs(end+1)=numcalculatedandcorezero-numcorezero*numcalculatedzero/91;
        allnumcorezero(end+1)=numcorezero;
        %sortedcalculatedfluxes=calculatedfluxes(sortedabscoreinds);
        
        %[corevalues calculatedfluxes sortedabscorevalues sortedabscoreinds calculatedfluxes(sortedabscoreinds)]
        allprctiles(:,i)=prctile(sortedabscorevalues,0:5:100);
    end
end
allnumcorezero(end+1)=mean(allnumcorezero);
celllinesarray2{end+1}='Average';
figure('Position',[300, 300, 1200, 500]);
a=bar(1:59,allnumcorezero);
set(gca,'XTickLabel',celllinesarray2);
set(gca,'XTick',1:59);
set(gca,'Ylim',[0 40]);
xticklabel_rotate;
    %saveas(a,strcat('eMOMACorrsummarize2/',strcat(stringstograph{i},'.png')));
allprctiles
for i=1:size(allprctiles,1)
    median(allprctiles(i,:))
end