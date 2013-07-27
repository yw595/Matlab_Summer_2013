function modelreturn = constrainexchange4(modelreturnoriginal,cellline)
[excnumarray exctextarray raw]=xlsread('../Supp Table 3 A community-driven global reconstruction of human metabolism 95.xls');
[height width]=size(excnumarray);
metsarray=exctextarray(10:149,2);
uniquemetstoexcrxninds=metstoexcrxns(metsarray,modelreturnoriginal,1);
celllinesarray=exctextarray(9,10:2:128);

exctokeep={'EX_gly(e)','EX_arg_L(e)','EX_asp_L(e)','EX_asn_L(e)','EX_cys_L(e)','EX_glu_L(e)','EX_gln_L(e)','EX_his_L(e)','EX_4hpro(e)','EX_ile_L(e)','EX_leu_L(e)',...
'EX_lys_L(e)','EX_met_L(e)','EX_phe_L(e)','EX_pro_L(e)','EX_ser_L(e)','EX_thr_L(e)','EX_trp_L(e)','EX_tyr_L(e)','EX_val_L(e)','EX_btn(e)','EX_chol(e)','EX_pnto_R(e)',...
'EX_fol(e)','EX_ncam(e)','EX_bz(e)','EX_pydxn(e)','EX_ribflv(e)','EX_thm(e)','EX_adpcbl(e)','EX_inost(e)','EX_ca2(e)','EX_so4(e)','EX_k(e)','EX_cl(e)','EX_na1(e)',...
'EX_hco3(e)','EX_pi(e)','EX_glc(e)','EX_gthrd(e)','EX_o2(e)','EX_co2(e)','EX_h2o(e)','EX_h(e)'};
celllineind=find(strcmp(cellline,celllinesarray));
modelreturn=modelreturnoriginal;
for i=1:length(modelreturn.rxns)
    if (~isempty(regexp(modelreturn.rxns{i},'^EX(.)*\(e\)$')))
        modelreturn.lb(i)=0;
    end
end
excnumarrayconstraintscol=excnumarray(:,celllineind*2+6);
for j=1:length(excnumarrayconstraintscol)
    met=metsarray{j};
    j
    excrxnind=uniquemetstoexcrxninds(met)
    if(excnumarrayconstraintscol(j)>0)
        for k=1:length(excrxnind)
            %modelreturn.ub(excrxnind(k))=modelreturn.ub(excrxnind(k))+excnumarrayconstraintscol(j)/length(excrxnind);
        end
    elseif(excnumarrayconstraintscol(j)<0)
        for k=1:length(excrxnind)
            if(sum(find(ismember(exctokeep,modelreturn.rxns{excrxnind(k)})))~=0)
                modelreturn.lb(excrxnind(k))=modelreturn.lb(excrxnind(k))+excnumarrayconstraintscol(j)/length(excrxnind);
            end
        end
    end
end
for j=1:length(exctokeep)
    ind=find(ismember(modelreturn.rxns,exctokeep{j}));
    if(modelreturn.lb(ind)==0)
        modelreturn.lb(ind)=-1000;
    end
    disp([num2str(ind) ' ' num2str(modelreturn.rxns{ind}) ' ' num2str(modelreturn.lb(ind)) ' ' num2str(modelreturn.ub(ind))]);
en
end

