%[excnumarray exctextarray raw]=xlsread('../Supp Table 3 A community-driven global reconstruction of human metabolism 95.xls');
%celllinesarray=exctextarray(9,10:2:128);
%for i=1:length(celllinesarray)
    %if(~strcmp(celllinesarray{i},'MDA-MB-468') && ~strcmp(celllinesarray{i},'RXF 393'))
constrainedfluxessofar={};
constrainedfluxesarray={};
glucosefluxesarray=[];
averagevaluesarray=[];
numfluxvaluesarray=[];
for i=1:92
    expressionFile=convertexpressionfilename('UACC-257');
    inputfile=['../eMOMACorroutconstrainedsequentialnarayanbetainescalejustFBS3/' expressionFile num2str(i) 'out'];
    %inputfile
    inputFI=fopen(inputfile,'r');

    line=fgetl(inputFI);
    inputvals=[];
    numnonzero=0;
    totalvalue=0;
    flagexfluxes=0;
    constrainedfluxes={};
    glucoseflux='';
    while line~=-1
        if(~isempty(regexp(line,'glucose')) && isempty(regexp(line,'UDP')))
            startindex=regexp(line,'(\-|\d|\.)+$');
            glucoseflux=line(startindex:end);
        end
        if(~isempty(regexp(line,'EX')))                				 
            %startindex=regexp(line,'(\-|\d|\.)+$');
	    startindex1=regexp(line,'(\-|\d|\.)+\t(\-|\d|\.)+$');
            startindex2=regexp(line,'(\-|\d|\.)+$');
	    excrxnname=line(1:startindex1-2);
            excrxnind=find(ismember(rec2.rxns,excrxnname));
            originallb=rec2.lb(excrxnind);
	    originalub=rec2.ub(excrxnind);
	    str2num(line(startindex1:startindex2-2));
	    if(str2num(line(startindex2:end))~=originalub)
	      %disp('HERE');
	    end
            if((str2num(line(startindex1:startindex2-2))~=originallb && str2num(line(startindex1:startindex2-2))~=0)|| str2num(line(startindex2:end))~=originalub)
                constrainedfluxes=union(constrainedfluxes,{line(1:startindex1-2)});
            end
        end
        if(flagexfluxes)
            startindex=regexp(line,'(\-|\d|\.)+$');
            exflux=str2num(line(startindex:length(line)));
            if(abs(exflux)>=.01)
                numnonzero=numnonzero+1;
                totalvalue=totalvalue+abs(exflux);
            end
        end
        if(~isempty(regexp(line,'v_solex')))
            flagexfluxes=1;
        end
        line=fgetl(inputFI);
    end
    %constrainedfluxes
    ismemberinds=ismember(constrainedfluxes,constrainedfluxessofar);
    if(isempty(constrainedfluxessofar) && ~isempty(constrainedfluxes))
        constrainedfluxesarray{end+1}='None';
        glucosefluxesarray(end+1)=-str2num(glucoseflux);
        averagevaluesarray(end+1)=totalvalue/numnonzero;
        numfluxvaluesarray(end+1)=numnonzero;
        disp([inputfile ' ' glucoseflux ' ' num2str(numnonzero) ' ' num2str(totalvalue/numnonzero)]);
    elseif(~isempty(constrainedfluxes(~ismemberinds)))
        temp=constrainedfluxes(~ismemberinds);
        newconstrainedflux=temp{1};
        constrainedfluxesarray{end+1}=newconstrainedflux;
        glucosefluxesarray(end+1)=-str2num(glucoseflux);
        averagevaluesarray(end+1)=totalvalue/numnonzero;
        numfluxvaluesarray(end+1)=numnonzero;
        disp([inputfile ' ' newconstrainedflux ' ' glucoseflux ' ' num2str(numnonzero) ' ' num2str(totalvalue/numnonzero)]);
    end
    %length(constrainedfluxessofar)
    constrainedfluxessofar=union(constrainedfluxessofar, constrainedfluxes);
    %if(~isempty(regexp(line,'glucose')) && isempty(regexp(line,'UDP')))                				 
    %disp([inputfile ' ' glucoseflux]);
    %end
    %disp([inputfile ' ' num2str(numnonzero) ' ' num2str(totalvalue/numnonzero)]);
end
figure('Position',[300, 300, 1200, 500],'Visible','off');
a=bar(1:length(constrainedfluxesarray),numfluxvaluesarray);
rubbish={};
for i=1:length(constrainedfluxesarray)
    rubbish{end+1}='a';
end
set(gca,'XTickLabel',constrainedfluxesarray);
set(gca,'XTick',1:length(constrainedfluxesarray));
set(gca,'Ylim',[0 30]);
title('Num fluxes');
xticklabel_rotate;
saveas(a,'../numfluxes.png');
	


