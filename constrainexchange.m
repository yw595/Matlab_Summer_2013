function model2 = constrainexchange(model1)
[excnumarray exctextarray raw]=xlsread('../Supp Table 3 A community-driven global reconstruction of human metabolism 95.xls');
[height width]=size(excnumarray);
jainmetsarray=exctextarray(10:100,1);
metsarray=exctextarray(10:100,2);
jainmetstomets=containers.Map(jainmetsarray,metsarray);
uniquemetsarray=values(jainmetstomets);
uniquemetstoexcrxnnames=metstoexcrxns(metsarray,model1,2);

exctokeep={'EX_gly(e)','EX_arg_L(e)','EX_asp_L(e)','EX_asn_L(e)','EX_cys_L(e)','EX_glu_L(e)','EX_gln_L(e)','EX_his_L(e)','EX_4hpro(e)','EX_ile_L(e)','EX_leu_L(e)',...
'EX_lys_L(e)','EX_met_L(e)','EX_phe_L(e)','EX_pro_L(e)','EX_ser_L(e)','EX_thr_L(e)','EX_trp_L(e)','EX_tyr_L(e)','EX_val_L(e)','EX_btn(e)','EX_chol(e)','EX_pnto_R(e)',...
'EX_fol(e)','EX_ncam(e)','EX_bz(e)','EX_pydxn(e)','EX_ribflv(e)','EX_thm(e)','EX_adpcbl(e)','EX_inost(e)','EX_ca2(e)','EX_so4(e)','EX_k(e)','EX_cl(e)','EX_na1(e)',...
'EX_hco3(e)','EX_pi(e)','EX_glc(e)','EX_gthrd(e)','EX_o2(e)','EX_co2(e)','EX_h2o(e)','EX_h(e)'};
model2=model1;
exctokeep2={};
for i=1:length(uniquemetsarray)
    excrxnnames=uniquemetstoexcrxnnames(uniquemetsarray{i});
    for j=1:length(excrxnnames)
        exctokeep2{end+1}=excrxnnames{j};
    end
end
for i=1:length(model2.rxns)
    if (~isempty(regexp(model2.rxns{i},'^EX(.)*\(e\)$')))
        if (sum(strcmp(model2.rxns{i},exctokeep))==0 && sum(strcmp(model2.rxns{i},exctokeep2))==0)
            model2.lb(i)=0;
        end
    end
end

