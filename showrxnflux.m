function showrxnflux(model,v_sol,rxns)
for i=1:length(rxns)
    rxn=rxns{i};
    chemicaleq=chemicalequation(model,rxn);
    rxnind=find(ismember(model.rxns,rxn));
    if(v_sol(rxnind)~=0)
        disp(sprintf([rxn '\t' model.rxnNames{rxnind} '\t' num2str(v_sol(rxnind)) '\t' chemicaleq]));
    end
end
end
