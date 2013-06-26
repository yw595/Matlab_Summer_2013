[excnumarray exctextarray raw]=xlsread('Supp Table 3 A community-driven global reconstruction of human metabolism 95.xls');
[height width]=size(excnumarray);
subexcnumarray=excnumarray(8:98,8:width);
jainmetsarray=exctextarray(10:100,1);
metsarray=exctextarray(10:100,2);
jainmetstomets=containers.Map(jainmetsarray,metsarray);
uniquemetsarray=values(jainmetstomets);
celllinesarray=exctextarray(9,10:2:128);
uniquemetstorxninds=containers.Map;
model=rec2;

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
        expressionFile=strrep(expressionarray{i},'(','_');
        expressionFile=strrep(expressionFile,')','_');
        expressionFile=strrep(expressionFile,' ','_');
        expressionFile=strrep(expressionFile,'/','_');
        expressionFile=strrep(expressionFile,'-','_');
        outputFile=strcat(strcat('eMOMACorrout2/',expressionFile),'out');
        expressionFile=strcat('NCI60exp/',strcat(expressionFile,'.csv'));
        
        runeMOMA(model,expressionFile,outputFile);
    end
end
