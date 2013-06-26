excarray=xlsread('Supp Table 3 A community-driven global reconstruction of human metabolism 95.xls');
[height width]=size(excarray);
subexcarray=excarray(8:98,8:width);
jainmetsarray={'carnosine','carnitine','homoserine','sucrose','thymidine','UDP-galactose/UDP-glucose'...
'thyroxine','triiodothyronine','ascorbate','DHAP','thiamine','biotin','glucuronate','GMP'...
'IMP','lithocholate','taurolithocholate','niacinamide','CMP','niacin','4-pyridoxate','folate'...
'tryptophan','anthranilate','serotonin','bilirubin','glycocholate','taurocholate',...
'taurodeoxycholate/taurochenodeoxycholate','choline','phosphocholine','betaine','spermine',...
'glutathione oxidized','lactose','isoleucine','5''-adenosylhomocysteine','methionine',...
'pantothenate','lysine','OH-phenylpyruvate','GABA','acetoacetate','adenine','adenosine',...
'alpha-ketoglutarate','alanine','AMP','arginine','asparagine','aspartate','citrate','citrulline',...
'creatine','creatinine','cytidine','2''-deoxyadenosine','2''-deoxycytidine','2''-deoxyuridine',...
'glucose','glutamine','glutamate','glycine','hypoxanthine','inosine','lactate','leucine','malate',...
'ornithine','orotate','oxalate','phenylalanine','proprionate','proline','pyruvate','sorbitol',...
'serine','spermidine','succinate','taurine','threonine','thymine','tyrosine','UMP','uracil',...
'urate','uridine','valine','citrate/isocitrate','glycerol_1','glycerol_2'};
metsarray={'carn','crn','hom_L','sucr','thymd','udpgal/udpg','thyox_L','triodthy',...
'ascb_L','dhap','thm','btn','glcur','gmp','imp','HC02191','HC02192','ncam','cmp',...
'nac','4pyrdx','fol','trp_L','anth','srtn','bilirub','gchola','tchola','tdchola',...
'chol','pchol_hs','glyb','sprm','gthox','lcts','ile_L','ahcys','met_L','pnto_R',...
'lys_L','34hpp','4abut','acac','ade','adn','akg','ala_L','amp','arg_L','asn_L',...
'asp_L','cit','citr_L','creat','creat','cytd','dad_2','dcyt','duri','glc_D','gln_L','glu_L',...
'gly','hxan','ins','lac_D/lac_L','leu_L','mal_L','orn','orot','oxa','phe_L','ppa',...
'pro_L','pyr','sbt_D','ser_L','spmd','succ','taur','thr_L','thym','tyr_L','ump','ura',...
'urate','uri','val_L','cit/icit','glyc','glyc'};
jainmetstomets=containers.Map(jainmetsarray,metsarray);
uniquemetsarray=values(jainmetstomets);
expressionarray={'UACC-257','OVCAR-8','OVCAR-5','SF-295','A549/ATCC','RXF 393',...
'CAKI-1','MDA-MB-435','EKVX','NCI-H23','HCC-2998','NCI_ADR-RES','HT29','M14',...
'HCT-116','HOP-62','786-O','NCI-H322M','BT-549','MALME-3M','TK-10','UO-31','T-47D',...
'MCF7','NCI-H460','SK-OV-3','SK-MEL-5','SNB-19','SF-539','COLO 205','KM12',...
'OVCAR-3','SF-268','DU-145','HS 578T','SW620','SK-MEL-28','A498','IGROV1','MOLT-4',...
'K562','HOP-92','PC-3','HL-60(TB)','MDA-MB-468','SK-MEL-2','NCI-H522','HCT-15','ACHN',...
'U251','LOX IMVI','SNB-75','OVCAR-4','NCI-H226','UACC-62','SN12C','MDA-MB-231_ATCC',...
'CCRF-CEM','RPMI 8226','SR'};
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
for i=1:length(expressionarray)
    if(~strcmp(expressionarray{i},'MDA-MB-468') && ~strcmp(expressionarray{i},'RXF 393'))
    expressionFile=strrep(expressionarray{i},'(','_');
    expressionFile=strrep(expressionFile,')','_');
    expressionFile=strrep(expressionFile,' ','_');
    expressionFile=strrep(expressionFile,'/','_');
    expressionFile=strrep(expressionFile,'-','_');
    outputfile=strcat(strcat('eMOMACorrout2/',expressionFile),'out');
    expressionFile=strcat('NCI60exp/',strcat(expressionFile,'.csv'));
    [modelIrrev,matchRev,rev2irrev,irrev2rev]=convertToIrreversible(model);
    [rxn_exp,rxn_exp_sd,rxn_rule_group]=computeMinDisj(modelIrrev,expressionFile);
    [v_solirrev, corrval, nvar]=eMOMA6Yiping(modelIrrev,rxn_exp,rxn_exp_sd,rxn_rule_group,{'pyruvate kinase'},0,{1});
    system(['touch ', outputfile]);
    outputFI=fopen(outputfile,'w');
    v_solex=zeros(length(jainmetsarray),1);
    fprintf(outputFI,'All fluxes from v_sol:\n');
    for j=1:length(v_solirrev)
        fprintf(outputFI,'%d\t%f\n',j,v_solirrev(j));
    end
    v_solrev=zeros(length(model.rxns),1);
    fprintf(outputFI,'All fluxes from v_solrev:\n');
    for j=1:length(irrev2rev)
        irrevrxnname=modelIrrev.rxns{j};
        if~isempty(regexp(irrevrxnname,'_b$'))
            v_solrev(irrev2rev(j))=v_solrev(irrev2rev(j))-v_solirrev(j);
        elseif~isempty(regexp(irrevrxnname,'_f$'))
            v_solrev(irrev2rev(j))=v_solrev(irrev2rev(j))+v_solirrev(j);
        elseif~isempty(regexp(irrevrxnname,'_r$'))
            v_solrev(irrev2rev(j))=-v_solirrev(j);
        else
            v_solrev(irrev2rev(j))=v_solirrev(j);
        end
    end
    for j=1:length(v_solrev)
        fprintf(outputFI,'%d\t%f\n',j,v_solrev(j));
    end
    fprintf(outputFI,'All fluxes from v_solex:\n');
    for j=1:length(v_solex)
        met=jainmetstomets(jainmetsarray{j});
        rxninds=uniquemetstorxninds(met);
        v_solex(j)=sum(v_solrev(rxninds));
        fprintf(outputFI,'%s\t%f\n',jainmetsarray{j},v_solex(j));
    end
    end
end
