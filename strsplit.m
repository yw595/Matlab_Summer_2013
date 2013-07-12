function words = strsplit(string, delimiter)
inds=strfind(string, delimiter);
words={};
for i=1:length(inds)
    if(i==1)
        words{end+1}=string(1:inds(i));
    elseif(i==length(inds))
        words{end+1}=string(inds(i)+length(delimiter):end);
    else
        words{end+1}=string(inds(i)+length(delimiter):inds(i+1)-1);
    end
end
end

