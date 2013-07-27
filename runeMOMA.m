function [v_solirrev v_solrev]=runeMOMA(model,expressionFile,varargin)
[modelIrrev,matchRev,rev2irrev,irrev2rev]=convertToIrreversible(model);
[rxn_exp,rxn_exp_sd,rxn_rule_group]=computeMinDisj(modelIrrev,expressionFile);
gene_to_scale={};
upt_const={};
%if(length(varargin)>0)
%    varargin{1}
%    varargin{2}
%end
if(nargin>2)
    varargincellarray=varargin{1};
    gene_to_scale=varargincellarray{1};
    varargincellarray{2}
    %str2num(varargincellarray{2})
    upt_const=varargincellarray{2};
else
    gene_to_scale={'pyruvate kinase'};
    upt_const={1};
end
[v_solirrev, corrval, nvar]=eMOMA6Yiping(modelIrrev,rxn_exp,rxn_exp_sd,rxn_rule_group,gene_to_scale,0,upt_const);

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