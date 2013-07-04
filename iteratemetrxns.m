function rxnstosearch = iteratemetrxns(model,v_sol)
metstosearch={'glc_d[c]','g6p[c]','f6p[c]','fdp[c]','dgap[c]','g3p[c]','13gdp[c]','3pg[c]','2pg[c]','pep[c]','pyr[c]',};
%metstosearch={'pep[c]','pyr[c]',};
rxnstosearch={};
for i=1:length(metstosearch)
    rxnstosearch=union(rxnstosearch,metrxns(model,v_sol,metstosearch{i}));
end