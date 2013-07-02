model=constrainexchange(rec2);
expressionFile='../A549_ATCC.csv2';
outputFile='../A549_ATCCout2';
outputFI=fopen(outputFile,'w');

[v_solirrev v_solrev]=runeMOMA(model,expressionFile);

fprintf(outputFI,'All fluxes from v_solirrev:\n');
for j=1:length(v_solirrev)
    fprintf(outputFI,'%d\t%f\n',j,v_solirrev(j));
end
fprintf(outputFI,'All fluxes from v_solrev:\n');
for j=1:length(v_solrev)
    fprintf(outputFI,'%d\t%f\n',j,v_solrev(j));
end
fclose(outputFI);