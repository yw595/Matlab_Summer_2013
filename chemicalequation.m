function equation = chemicalequation(model,rxn)
rxnind=find(ismember(model.rxns,rxn));
rxnmetsinds=find(model.S(:,rxnind));
rxnmetscoefs=model.S(rxnmetsinds,rxnind);
reactants='';
products='';
for i=1:length(rxnmetsinds)
    if(rxnmetscoefs(i)<0)
        reactants=[reactants, num2str(rxnmetscoefs(i)*-1), ' ', model.mets{rxnmetsinds(i)}, ' '];
    else
        products=[products num2str(rxnmetscoefs(i)) ' ' model.mets{rxnmetsinds(i)} ' '];
    end
end
equation=[strtrim(reactants) ' => ' strtrim(products)];
end

