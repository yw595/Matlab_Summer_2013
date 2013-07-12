[excnumarray exctextarray raw]=xlsread('../Supp Table 3 A community-driven global reconstruction of human metabolism 95.xls');
[height width]=size(excnumarray);
jainmetsarray=exctextarray(10:100,1);
metsarray=exctextarray(10:100,2);
jainmetstomets=containers.Map(jainmetsarray,metsarray);
model=constrainexchange(rec2);
uniquemetstorxninds=metstoexcrxns(metsarray,model,1);
celllinesarray=exctextarray(9,10:2:128);

for i=1:length(celllinesarray)
    if(~strcmp(celllinesarray{i},'MDA-MB-468') && ~strcmp(celllinesarray{i},'RXF 393'))
        expressionFile=convertexpressionfilename(celllinesarray{i});
        outputFile=['../eMOMACorroutconstraineddummy/' expressionFile 'out'];
        expressionFileLoc=['../NCI60exp/' expressionFile '.csv'];
        outputFI=fopen(outputFile,'w');
        
        specificmodel=constrainexchange2(model,celllinesarray{i});
        [v_solirrev v_solrev]=runeMOMA(specificmodel,expressionFileLoc);
        
        fprintf(outputFI,'All fluxes from v_solirrev:\n');
        for j=1:length(v_solirrev)
            fprintf(outputFI,'%d\t%f\n',j,v_solirrev(j));
        end
        fprintf(outputFI,'All fluxes from v_solrev:\n');
        for j=1:length(v_solrev)
            fprintf(outputFI,'%d\t%f\n',j,v_solrev(j));
        end
        v_solex=zeros(length(jainmetsarray),1);
        fprintf(outputFI,'All fluxes from v_solex:\n');
        for j=1:length(v_solex)
            met=jainmetstomets(jainmetsarray{j});
            rxninds=uniquemetstorxninds(met);
            v_solex(j)=sum(v_solrev(rxninds));
            disp(sprintf('%s\t%f',jainmetsarray{j},v_solex(j)));
            fprintf(outputFI,'%s\t%f\n',jainmetsarray{j},v_solex(j));
        end
        fclose(outputFI);
    end
end