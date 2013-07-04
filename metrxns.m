function rxnsreturn=metrxns(model,v_sol,met)
metind=find(ismember(model.mets,met));
metrxninds=find(model.S(metind,:));
rxnsreturn={};
for i=1:length(model.rxns(metrxninds))
    rxn=model.rxns{metrxninds(i)};
    %chemicaleq=chemicalequation(model,rxn);
    %rxnind=find(ismember(model.rxns,rxn));
    if(v_sol(metrxninds(i))~=0)
        rxnsreturn{end+1}=rxn;
    end
end
end

