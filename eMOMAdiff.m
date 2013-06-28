function [v_soldiffabssorted sortInds] = eMOMAdiff(inputFile1, inputFile2)
inputFI1=fopen(inputFile1,'r');
inputFI2=fopen(inputFile2,'r');

line1=fgetl(inputFI1);
line2=fgetl(inputFI2);
inputvals=[];
flagrevfluxes=0;
v_soldiff=[];
while line1~=-1
    if(flagrevfluxes)
        startindex1=regexp(line1,'(\-|\d|\.)+$');
        startindex2=regexp(line2,'(\-|\d|\.)+$');
        v_soldiff(end+1)=str2num(line1(startindex1:length(line1)))-str2num(line2(startindex2:length(line2)));
    end
    if(~isempty(regexp(line1,'v_solrev')))
        flagrevfluxes=1;
    end
    line1=fgetl(inputFI1);
    line2=fgetl(inputFI2);
end
[v_soldiffabssorted sortInds]=sort(abs(v_soldiff));
end

