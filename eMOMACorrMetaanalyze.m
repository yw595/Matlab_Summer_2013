expressionarray={'UACC-257','OVCAR-8','OVCAR-5','SF-295','A549/ATCC','RXF 393',...
'CAKI-1','MDA-MB-435','EKVX','NCI-H23','HCC-2998','NCI_ADR-RES','HT29','M14',...
'HCT-116','HOP-62','786-O','NCI-H322M','BT-549','MALME-3M','TK-10','UO-31','T-47D',...
'MCF7','NCI-H460','SK-OV-3','SK-MEL-5','SNB-19','SF-539','COLO 205','KM12',...
'OVCAR-3','SF-268','DU-145','HS 578T','SW620','SK-MEL-28','A498','IGROV1','MOLT-4',...
'K562','HOP-92','PC-3','HL-60(TB)','MDA-MB-468','SK-MEL-2','NCI-H522','HCT-15','ACHN',...
'U251','LOX IMVI','SNB-75','OVCAR-4','NCI-H226','UACC-62','SN12C','MDA-MB-231_ATCC',...
'CCRF-CEM','RPMI 8226','SR'};
expressionarray2={};
for i=1:length(expressionarray)
    if(~strcmp(expressionarray{i},'MDA-MB-468') && ~strcmp(expressionarray{i},'RXF 393'))
        expressionarray2{end+1}=expressionarray{i};
    end
end
input1='';
input2='';
inputFI1=fopen(input1,'r');
inputFI2=fopen(input2,'r');
pearsondiffs=[];
spearmandiffs=[];
kendalldiffs=[];
cossimdiffs=[];
avgfluxdiffs=[];
uptakesensdiffs=[];
releasesensdiffs=[];

line1=fgetl(inputFI1);
inputvals=[];
while line1~=-1
    if(~isempty(regexp(line1,'pearson')))
        startindex1=regexp(line1,'(\-|\d|\.)+$');
        startindex2=regexp(line2,'(\-|\d|\.)+$');
        diff=str2num(line1(startindex1:end))-str2num(line2(startindex2:end));
        pearsondiffs(end+1)=diff;
    elseif(~isempty(regexp(line1,'spearman')))
        startindex1=regexp(line1,'(\-|\d|\.)+$');
        startindex2=regexp(line2,'(\-|\d|\.)+$');
        diff=str2num(line1(startindex1:end))-str2num(line2(startindex2:end));
        spearmandiffs(end+1)=diff;
    elseif(~isempty(regexp(line1,'kendall')))
        startindex1=regexp(line1,'(\-|\d|\.)+$');
        startindex2=regexp(line2,'(\-|\d|\.)+$');
        diff=str2num(line1(startindex1:end))-str2num(line2(startindex2:end));
        kendalldiffs(end+1)=diff;
    elseif(~isempty(regexp(line1,'cosine')))
        startindex1=regexp(line1,'(\-|\d|\.)+$');
        startindex2=regexp(line2,'(\-|\d|\.)+$');
        diff=str2num(line1(startindex1:end))-str2num(line2(startindex2:end));
        cossimdiffs(end+1)=diff;
    elseif(~isempty(regexp(line1,'avg flux')))
        startindex1=regexp(line1,'(\-|\d|\.)+$');
        startindex2=regexp(line2,'(\-|\d|\.)+$');
        diff=str2num(line1(startindex1:end))-str2num(line2(startindex2:end));
        avgfluxdiffdiffs(end+1)=diff;
    elseif(~isempty(regexp(line1,'uptakesens')))
        startindex1=regexp(line1,'(\-|\d|\.)+$');
        startindex2=regexp(line2,'(\-|\d|\.)+$');
        diff=str2num(line1(startindex1:end))-str2num(line2(startindex2:end));
        uptakesensdiffs(end+1)=diff;
    elseif(~isempty(regexp(line1,'releasesens')))
        startindex1=regexp(line1,'(\-|\d|\.)+$');
        startindex2=regexp(line2,'(\-|\d|\.)+$');
        diff=str2num(line1(startindex1:end))-str2num(line2(startindex2:end));
        releasesensdiffs(end+1)=diff;
    end
    line1=fgetl(inputFI1);
    line2=fgetl(inputFI2);
end
pearsondiffs(end+1)=mean(pearsondiffs);
spearmandiffs(end+1)=mean(spearmandiffs);
kendalldiffs(end+1)=mean(kendalldiffs);
cossimdiffs(end+1)=mean(cossimdiffs);
avgfluxdiffdiffs(end+1)=mean(avgfluxdiffdiffs);
uptakesensdiffs(end+1)=mean(releasesensdiffs);
releasesensdiffs(end+1)=mean(uptakesensdiffs);
expressionarray2{end+1}='Average';