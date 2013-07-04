[excnumarray exctextarray raw]=xlsread('../Supp Table 3 A community-driven global reconstruction of human metabolism 95.xls');
[height width]=size(excnumarray);
subexcnumarray=excnumarray(8:98,8:width);
jainmetsarray=exctextarray(10:100,1);
metsarray=exctextarray(10:100,2);
jainmetstomets=containers.Map(jainmetsarray,metsarray);
uniquemetsarray=values(jainmetstomets);
uniquemetstorxninds=containers.Map;
celllinesarray=exctextarray(9,10:2:128);
model=constrainexchange(rec2);

for i=1:length(metsarray)
    if(strcmp(metsarray{i},'34hpp'))
        excrxnname='EX_34hpp';
        excrxnind=find(strcmp(excrxnname,model.rxns));
        uniquemetstorxninds(metsarray{i})=excrxnind;
    elseif(strcmp(metsarray{i},'glc_D'))
        excrxnname='EX_glc(e)';
        excrxnind=find(strcmp(excrxnname,model.rxns));
        uniquemetstorxninds(metsarray{i})=excrxnind;
    elseif(strcmp(metsarray{i},'udpgal/udpg'))
        excrxnname1='EX_udpgal(e)';
        excrxnind1=find(strcmp(excrxnname1,model.rxns));
        excrxnname2='EX_udpg(e)';
        excrxnind2=find(strcmp(excrxnname2,model.rxns));
        uniquemetstorxninds(metsarray{i})=[excrxnind1(1) excrxnind2(1)];
    elseif(strcmp(metsarray{i},'lac_D/lac_L'))
        excrxnname1='EX_lac_D(e)';
        excrxnind1=find(strcmp(excrxnname1,model.rxns));
        excrxnname2='EX_lac_L(e)';
        excrxnind2=find(strcmp(excrxnname2,model.rxns));
        uniquemetstorxninds(metsarray{i})=[excrxnind1(1) excrxnind2(1)];
    elseif(strcmp(metsarray{i},'sbt_D'))
        excrxnname='EX_sbt-d(e)';
        excrxnind=find(strcmp(excrxnname,model.rxns));
        uniquemetstorxninds(metsarray{i})=excrxnind;
    else
        excrxnname=strcat('EX_',strcat(metsarray{i},'(e)'));
        excrxnind=find(strcmp(excrxnname,model.rxns));
        uniquemetstorxninds(metsarray{i})=excrxnind;
    end
end

for i=1:length(celllinesarray)
    if(~strcmp(celllinesarray{i},'MDA-MB-468') && ~strcmp(celllinesarray{i},'RXF 393'))
        expressionFile=strrep(celllinesarray{i},'(','_');
        expressionFile=strrep(expressionFile,')','_');
        expressionFile=strrep(expressionFile,' ','_');
        expressionFile=strrep(expressionFile,'/','_');
        expressionFile=strrep(expressionFile,'-','_');
        outputFile=['../eMOMACorroutconstrainedextra/' expressionFile 'out'];
        expressionFile=['../NCI60exp/' expressionFile '.csv'];
        outputFI=fopen(outputFile,'w');
        
        [v_solirrev v_solrev]=runeMOMA(model,expressionFile);
        
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
            fprintf(outputFI,'%s\t%f\n',jainmetsarray{j},v_solex(j));
        end
        fclose(outputFI);
    end
end