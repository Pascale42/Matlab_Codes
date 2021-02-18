function Par = LoadParCSD(FileName) 

Arr = LoadStringArray(FileName);
Par.nShanks = size(Arr,1);
Par.Channels=cell(Par.nShanks,1);
for i=1:Par.nShanks
    Par.nChannels(i) = 0;
    Par.Location = Arr{i,1};
    for j=2:size(Arr,2)
        if ~isempty(Arr{i,j})
            Par.Channels{i}(end+1) = str2num(Arr{i,j})+1; % adds 1 so now channels are from 1!!
            Par.nChannels(i) =     Par.nChannels(i) +1;
        end
    end
end
