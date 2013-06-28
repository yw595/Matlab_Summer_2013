function [v_solirrev v_solrev]=runeMOMA(model,expressionFile)
[modelIrrev,matchRev,rev2irrev,irrev2rev]=convertToIrreversible(model);
[rxn_exp,rxn_exp_sd,rxn_rule_group]=computeMinDisj(modelIrrev,expressionFile);
[v_solirrev, corrval, nvar]=eMOMA6Yiping(modelIrrev,rxn_exp,rxn_exp_sd,rxn_rule_group,{'pyruvate kinase'},0,{1});

v_solrev=zeros(length(model.rxns),1);
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
end