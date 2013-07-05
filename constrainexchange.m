function model2 = constrainexchange(model1)
[excnumarray exctextarray raw]=xlsread('../Supp Table 3 A community-driven global reconstruction of human metabolism 95.xls');
[height width]=size(excnumarray);
jainmetsarray=exctextarray(10:100,1);
metsarray=exctextarray(10:100,2);
jainmetstomets=containers.Map(jainmetsarray,metsarray);
uniquemetsarray=values(jainmetstomets);
uniquemetstoexcrxnnames=containers.Map;

for i=1:length(uniquemetsarray)
    if(strcmp(uniquemetsarray{i},'34hpp'))
        excrxnname='EX_34hpp';
        uniquemetstoexcrxnnames(uniquemetsarray{i})={excrxnname};
    elseif(strcmp(uniquemetsarray{i},'glc_D'))
        excrxnname='EX_glc(e)';
        uniquemetstoexcrxnnames(uniquemetsarray{i})={excrxnname};
    elseif(strcmp(uniquemetsarray{i},'udpgal/udpg'))
        excrxnname1='EX_udpgal(e)';
        excrxnname2='EX_udpg(e)';
        uniquemetstoexcrxnnames(uniquemetsarray{i})={excrxnname1,excrxnname2};
        %uniquemetstoexcrxnnames(uniquemetsarray{i})
    elseif(strcmp(uniquemetsarray{i},'lac_D/lac_L'))
        excrxnname1='EX_lac_D(e)';
        excrxnname2='EX_lac_L(e)';
        uniquemetstoexcrxnnames(uniquemetsarray{i})={excrxnname1,excrxnname2};
    elseif(strcmp(uniquemetsarray{i},'sbt_D'))
        excrxnname='EX_sbt-d(e)';
        uniquemetstoexcrxnnames(uniquemetsarray{i})={excrxnname};
    elseif(strcmp(uniquemetsarray{i},'tyr_l'))
        excrxnname='EX_tyr_L(e)';
        uniquemetstoexcrxnnames(uniquemetsarray{i})={excrxnname};
    elseif(strcmp(uniquemetsarray{i},'cit/icit'))
        excrxnname='EX_cit(e)';
        uniquemetstoexcrxnnames(uniquemetsarray{i})={excrxnname};
    else
        excrxnname=strcat('EX_',strcat(uniquemetsarray{i},'(e)'));
        uniquemetstoexcrxnnames(uniquemetsarray{i})={excrxnname};
    end
end
exctokeep={'EX_gly(e)','EX_arg_L(e)','EX_asp_L(e)','EX_asn_L(e)','EX_cys_L(e)','EX_glu_L(e)','EX_gln_L(e)','EX_his_L(e)','EX_4hpro(e)','EX_ile_L(e)','EX_leu_L(e)',...
'EX_lys_L(e)','EX_met_L(e)','EX_phe_L(e)','EX_pro_L(e)','EX_ser_L(e)','EX_thr_L(e)','EX_trp_L(e)','EX_tyr_L(e)','EX_val_L(e)','EX_btn(e)','EX_chol(e)','EX_pnto_R(e)',...
'EX_fol(e)','EX_ncam(e)','EX_bz(e)','EX_pydxn(e)','EX_ribflv(e)','EX_thm(e)','EX_adpcbl(e)','EX_inost(e)','EX_ca2(e)','EX_so4(e)','EX_k(e)','EX_cl(e)','EX_na1(e)',...
'EX_hco3(e)','EX_pi(e)','EX_glc(e)','EX_gthrd(e)','EX_o2(e)','EX_co2(e)','EX_h2o(e)','EX_h(e)'};
%exctokeep2={'EX_carn(e)','EX_crn(e)','EX_hom_L(e)','EX_sucr(e)','EX_thymd(e)','EX_udpgal(e)','EX_udpg(e)','EX_thyox_L(e)','EX_triodthy(e)','EX_ascb_L(e)','EX_dhap(e)',...
%'EX_thm(e)','EX_btn(e)','EX_glcur(e)','EX_gmp(e)','EX_imp(e)','EX_HC02191(e)','EX_HC02192(e)','EX_ncam(e)','EX_cmp(e)','EX_nac(e)','EX_pyrdx(e)',...
%'EX_anth(e)','EX_srtn(e)','EX_bilirub(e)','EX_gchola(e)','EX_tchola(e)','EX_tdchola(e)','EX_chol(e)','EX_pchol_hs(e)','EX_glyb(e)','EX_sprm(e)','EX_gthox(e)','EX_lcts(e)',...
%'EX_ancys(e)','EX_34hpp(e)','EX_4abut(e)','EX_acac(e)','EX_ade(e)','EX_adn(e)','EX_akg(e)','EX_amp(e)','EX_cit(e)','EX_citr_L(e)','EX_creat(e)','EX_cytd(e)',...
%'EX_dad_2(e)','EX_dcyt(e)','EX_duri(e)','EX_hxan(e)','EX_ins(e)','EX_lac_L(e)','EX_lac_D(e)','EX_mal_L(e)','EX_orn(e)','EX_orot(e)','EX_oxa(e)',...
%'EX_ppa(e)','EX_pyr(e)','EX_sbt_D(e)','EX_spmd(e)','EX_succ(e)','EX_taur(e)','EX_thym(e)','EX_ump(e)','EX_ura(e)','EX_urate(e)','EX_uri(e)','EX_cit(e)','EX_glyc(e)'};
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
for i=1:length(model2.rxns)
    if (~isempty(regexp(model2.rxns{i},'^EX(.)*\(e\)$')))
        if (model2.lb(i)~=0)
            model2.rxns{i}
        end
    end
end

