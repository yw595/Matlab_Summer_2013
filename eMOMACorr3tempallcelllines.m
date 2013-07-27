exctokeep={'EX_gly(e)','EX_arg_L(e)','EX_asp_L(e)','EX_asn_L(e)','EX_cys_L(e)','EX_glu_L(e)','EX_gln_L(e)','EX_his_L(e)','EX_4hpro(e)','EX_ile_L(e)','EX_leu_L(e)',...
'EX_lys_L(e)','EX_met_L(e)','EX_phe_L(e)','EX_pro_L(e)','EX_ser_L(e)','EX_thr_L(e)','EX_trp_L(e)','EX_tyr_L(e)','EX_val_L(e)','EX_btn(e)','EX_chol(e)','EX_pnto_R(e)',...
'EX_fol(e)','EX_ncam(e)','EX_bz(e)','EX_pydxn(e)','EX_ribflv(e)','EX_thm(e)','EX_adpcbl(e)','EX_inost(e)','EX_ca2(e)','EX_so4(e)','EX_k(e)','EX_cl(e)','EX_na1(e)',...
'EX_hco3(e)','EX_pi(e)','EX_glc(e)','EX_gthrd(e)','EX_o2(e)','EX_co2(e)','EX_h2o(e)','EX_h(e)'};
[excnumarray exctextarray raw]=xlsread('../Supp Table 3 A community-driven global reconstruction of human metabolism 95.xls');
[height width]=size(excnumarray);
jainmetsarray=exctextarray(10:100,1);
metsarray=exctextarray(10:100,2);
excrxnnames1=values(metstoexcrxns(metsarray,model,0));
excrxnnames={};
for i=1:length(excrxnnames1)
  temp=excrxnnames1{i};
  for j=1:length(temp)
    excrxnnames{end+1}=temp{j};
  end
end
excrxnnames;
jainmetstomets=containers.Map(jainmetsarray,metsarray);
model=constrainexchange(rec2);
uniquemetstorxninds=metstoexcrxns(metsarray,model,1);
celllinesarray=exctextarray(9,10:2:128);

for i=1:length(celllinesarray)
    if(~strcmp(celllinesarray{i},'MDA-MB-468') && ~strcmp(celllinesarray{i},'RXF 393'))
    models=constrainexchange3(rec2,celllinesarray{i});
    tempmodel=models{end};
    expressionFile='UACC_257';
    outputFile=['../eMOMACorroutconstrainedsequentialnarayanbetainescalejustFBS2/' expressionFile num2str(i) 'out'];
    expressionFileLoc=['../NCI60exp/' expressionFile '.csv'];
    outputFI=fopen(outputFile,'w');

    %if(i==1)
    %    [v_solirrev v_solrev]=runeMOMA(models{i},expressionFileLoc);
    %else
        %ithmodel=models{i};
        glcexcind=find(ismember(tempmodel.rxns,'EX_glyb(e)'));
        glcexclb=tempmodel.lb(glcexcind);
        glcexclb=glcexclb(1);
        [v_solirrev v_solrev]=runeMOMA(models{end},expressionFileLoc,{{'pyruvate kinase'},{abs(glcexclb)}});
    %end
        
	fprintf(outputFI,'All lower and upper bounds:\n');
	%tempmodel=models{i};
	for j=1:length(tempmodel.rxns)
        tempmodelrxns=tempmodel.rxns;
	%rxnind=find(ismember(tempmodelrxns{j},
        if(sum(strcmp(tempmodelrxns{j},exctokeep))~=0 || sum(strcmp(tempmodelrxns{j},excrxnnames))~=0)
            fprintf(outputFI,'%s\t%f\t%f\n',tempmodelrxns{j},tempmodel.lb(j),tempmodel.ub(j));
        end
	end
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