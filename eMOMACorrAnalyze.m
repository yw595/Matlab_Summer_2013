[exctextarray excnumarray raw]=xlsread('Supp Table 3 A community-driven global reconstruction of human metabolism 95.xls');
[height width]=size(excnumarray);
subexcnumarray=excnumarray(8:98,8:width);
celllinesarray=exctextarray(9,10:2:128);
outputfile='eMOMACorrsummarize2/summary.txt';
outputFI=fopen(outputfile,'w');
expressionarray2={};
allpearsonrho=[];
allspearmanrho=[];
allkendallrho=[];
allcossim=[];
allavgfluxdiff=[];
allreleasesens=[];
alluptakesens=[];
allnumreleasing=[];
allnumuptaking=[];
allnumzero=[];
allnumexcluded=[];
allnumexcludedfromfva=[];

for i=1:length(expressionarray)
    if(~strcmp(expressionarray{i},'MDA-MB-468') && ~strcmp(expressionarray{i},'RXF 393'))
        expressionFile=strrep(expressionarray{i},'(','_');
        expressionFile=strrep(expressionFile,')','_');
        expressionFile=strrep(expressionFile,' ','_');
        expressionFile=strrep(expressionFile,'/','_');
        expressionFile=strrep(expressionFile,'-','_');
        inputfile=strcat(strcat('eMOMACorrout2/',expressionFile),'out');
        inputFI=fopen(inputfile,'r');
        
        line=fgetl(inputFI);
        inputvals=[];
        flagexfluxes=0;
        while line~=-1
            if(flagexfluxes)
                startindex=regexp(line,'(\-|\d|\.)+$');
                inputvals=[inputvals; str2num(line(startindex:length(line)))];
            end
            if(~isempty(regexp(line,'v_solex')))
                flagexfluxes=1;
            end
            line=fgetl(inputFI);
        end
        
        subexcarrayabs=abs(subexcarray(:,i*2));
        subexcarraysign=sign(subexcarray(:,i*2));
        sortsubexcarrayabs=sort(subexcarrayabs);
        
        uptaketruepos=0;
        uptakefalseneg=0;
        releasetruepos=0;
        releasefalseneg=0;
        numreleasing=0;
        numuptaking=0;
        numzero=0;
        numexcluded=0;
        numexcludedfromfva=0;
        includedinds=[];
        
        for j=1:size(subexcarray(:,i*2),1)
            %if(abs(inputvals(j))<.00001)
            %    numexcluded=numexcluded+1;
            %if(subexcarray(j,i*2)>0 && excarray(j+7,3)==0)
            %    numexcludedfromfva=numexcludedfromfva+1;
            %elseif(subexcarray(j,i*2)<0 && excarray(j+7,1)==0)
            %    numexcludedfromfva=numexcludedfromfva+1;
            if(subexcarray(j,i*2)>0)
                numreleasing=numreleasing+1;
                includedinds(end+1)=j;
                if(inputvals(j)>0)
                    releasetruepos=releasetruepos+1;
                else
                    releasefalseneg=releasefalseneg+1;
                end
            elseif(subexcarray(j,i*2)<0)
                numuptaking=numuptaking+1;
                includedinds(end+1)=j;
                if(inputvals(j)<0)
                    uptaketruepos=uptaketruepos+1;
                else
                    uptakefalseneg=uptakefalseneg+1;
                end
            else
                numzero=numzero+1;
                includedinds(end+1)=j;
            end
        end
        %includedinds
        uptakefalsepos=releasefalseneg;
        uptaketrueneg=releasetruepos;
        releasefalsepos=uptakefalseneg;
        releasetrueneg=uptaketruepos;
        
        uptakesens=uptaketruepos/(uptaketruepos+uptakefalseneg);
        releasesens=releasetruepos/(releasetruepos+releasefalseneg);
        
        [emomapearsonrho emomapearsonpval]=corr(inputvals(includedinds),subexcarray(includedinds,i*2),'type','Pearson');
        [emomaspearmanrho emomaspearmanpval]=corr(inputvals(includedinds),subexcarray(includedinds,i*2),'type','Spearman');
        [emomakendallrho emomakendallpval]=corr(inputvals(includedinds),subexcarray(includedinds,i*2),'type','Kendall');
        cossim=(inputvals(includedinds)'*subexcarray(includedinds,i*2))/(norm(inputvals(includedinds))*norm(subexcarray(includedinds,i*2)));
        avgfluxdiff=mean(inputvals(includedinds)-subexcarray(includedinds,i*2));
        
        allpearsonrho(end+1)=emomapearsonrho;
        allspearmanrho(end+1)=emomaspearmanrho;
        allkendallrho(end+1)=emomakendallrho;
        allcossim(end+1)=cossim;
        allavgfluxdiff(end+1)=avgfluxdiff;
        allreleasesens(end+1)=releasesens;
        alluptakesens(end+1)=uptakesens;
        allnumreleasing(end+1)=numreleasing;
        allnumuptaking(end+1)=numuptaking;
        allnumzero(end+1)=numzero;
        allnumexcluded(end+1)=numexcluded;
        allnumexcludedfromfva(end+1)=numexcludedfromfva;
        expressionarray2{end+1}=expressionarray{i};
        
        fprintf(outputFI,'%s\n',expressionarray{i}); 
        fprintf(outputFI,'pearson corr: %f\n', emomapearsonrho);
        fprintf(outputFI,'spearman corr: %f\n', emomaspearmanrho);
        fprintf(outputFI,'kendall corr: %f\n', emomakendallrho);
        fprintf(outputFI,'cosine sim: %f\n', cossim);
        fprintf(outputFI,'avg flux diff: %f\n', avgfluxdiff);
        fprintf(outputFI,'uptakesens: %f\n', uptakesens);
        fprintf(outputFI,'releasesens: %f\n', releasesens);
        %fprintf(outputFI,'releasetruepos: %d releasefalseneg: %d uptaketruepos: %d uptakefalseneg: %d\n',releasetruepos,releasefalseneg,uptaketruepos,uptakefalseneg);
        fprintf(outputFI,'numreleasing: %d numuptaking: %d numzero: %d numexcluded: %d numexcludedfromfva: %d\n',numreleasing,numuptaking,numzero, numexcluded, numexcludedfromfva);
    end
end
allpearsonrho(end+1)=mean(allpearsonrho);
allspearmanrho(end+1)=mean(allspearmanrho);
allkendallrho(end+1)=mean(allkendallrho);
allcossim(end+1)=mean(allcossim);
allavgfluxdiff(end+1)=mean(allavgfluxdiff);
allreleasesens(end+1)=mean(allreleasesens);
alluptakesens(end+1)=mean(alluptakesens);
allnumreleasing(end+1)=mean(allnumreleasing);
allnumuptaking(end+1)=mean(allnumuptaking);
allnumzero(end+1)=mean(allnumzero);
allnumexcluded(end+1)=mean(allnumexcluded);
allnumexcludedfromfva(end+1)=mean(allnumexcludedfromfva);
expressionarray2{end+1}='Average';
fprintf(outputFI,'%s\n',expressionarray{end}); 
fprintf(outputFI,'pearson corr: %f\n', allpearsonrho{end});
fprintf(outputFI,'spearman corr: %f\n', allspearmanrho{end});
fprintf(outputFI,'kendall corr: %f\n', allkendallrho{end});
fprintf(outputFI,'cosine sim: %f\n', allcossim{end});
fprintf(outputFI,'avg flux diff: %f\n', allavgfluxdiff{end});
fprintf(outputFI,'uptakesens: %f\n', alluptakesens{end});
fprintf(outputFI,'releasesens: %f\n', allreleasesens{end});
%fprintf(outputFI,'releasetruepos: %d releasefalseneg: %d uptaketruepos: %d uptakefalseneg: %d\n',releasetruepos,releasefalseneg,uptaketruepos,uptakefalseneg);
fprintf(outputFI,'numreleasing: %d numuptaking: %d numzero: %d numexcluded: %d numexcludedfromfva: %d\n',...
allnumreleasing(end),allnumuptaking(end),allnumzero(end),allnumexcluded(end),allnumexcludedfromfva(end));

outputstograph={allpearsonrho,allspearmanrho,allkendallrho,allcossim,allavgfluxdiff,allreleasesens,alluptakesens,allnumreleasing,allnumzero};
stringstograph={'allpearsonrho','allspearmanrho','allkendallrho','allcossim','allavgfluxdiff','allreleasesens','alluptakesens','allnumreleasing','allnumzero'};
ylimstograph={[-1 1],[-1 1],[-1 1],[-1 1],[0 10],[0 1],[0 1],[0 50],[0 10]};
for i=1:length(outputstograph)
    figure('Position',[300, 300, 1200, 500],'Visible','off');
    a=bar(1:59,outputstograph{i});
    set(gca,'XTickLabel',expressionarray2);
    set(gca,'XTick',1:59);
    set(gca,'Ylim',ylimstograph{i});
    title(stringstograph{i});
    xticklabel_rotate;
    saveas(a,strcat('eMOMACorrsummarize2/',strcat(stringstograph{i},'.png')));
end



