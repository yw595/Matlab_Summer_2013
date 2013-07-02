function [v_solfinratiossorted v_solfindiffssorted finInds v_solinfdiffssorted infInds v_soldiffs v_solratios] = eMOMAdiff(inputFile1, inputFile2)
inputFI1=fopen(inputFile1,'r');
inputFI2=fopen(inputFile2,'r');

line1=fgetl(inputFI1);
line2=fgetl(inputFI2);
inputvals=[];
flagrevfluxes=0;
v_soldiffs=[];
v_solratios=[];
v_sol1=[];
v_sol2=[];
while line1~=-1
    if(flagrevfluxes)
        startindex1=regexp(line1,'(\-|\d|\.)+$');
        startindex2=regexp(line2,'(\-|\d|\.)+$');
        v_sol1(end+1)=str2num(line1(startindex1:length(line1)));
        v_sol2(end+1)=str2num(line2(startindex2:length(line2)));
        v_soldiffs(end+1)=str2num(line1(startindex1:length(line1)))-str2num(line2(startindex2:length(line2)));
        v_solratios(end+1)=(str2num(line1(startindex1:length(line1)))-str2num(line2(startindex2:length(line2))))/(str2num(line1(startindex1:length(line1))));
    end
    if(~isempty(regexp(line1,'v_solrev')))
        flagrevfluxes=1;
    end
    line1=fgetl(inputFI1);
    line2=fgetl(inputFI2);
end
v_solInds=1:7440;
infInds=v_solInds(isinf(v_solratios));
[v_solinfdiffssorted infIndssortInds]=sort(abs(v_soldiffs(infInds)));
infInds=infInds(infIndssortInds);
v_solinfratiossorted=v_solratios(infInds);

finInds=v_solInds(isfinite(v_solratios));
[v_solfinratiossorted finIndssortInds]=sort(abs(v_solratios(finInds)));
finInds=finInds(finIndssortInds);
v_solfinratiossorted=v_solratios(finInds);
v_solfindiffssorted=v_soldiffs(finInds);

%[v_solratiosabssorted sortInds]=sort(abs(v_solratios));
%v_solratiossorted=v_solratios(sortInds);
%finInds=isfinite(v_solratiossorted);
%v_solratiossorted=v_solratiossorted(finInds);
%sortInds=sortInds(finInds);
end

