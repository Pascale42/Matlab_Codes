%function FileCSD(FileBase, nChannels, Channels2Use, SpatialSmoother)
% does a CSD on a file, smoother with [1 2 1] by default
function FileCSD(FileBase, varargin)
FileIn = [FileBase '.eeg'];
FileOut = [FileBase '.csd'];
if FileExists([FileBase '.xml'])
    Par = LoadPar([FileBase '.xml']);
    nChannels = Par.nChannels;
    nShanks = length(Par.AnatGrps);
    for i=1:nShanks
        if length(Par.AnatGrps(i).Channels)<3; continue; end
        Channels2Use{i} = Par.AnatGrps(i).Channels(find(~Par.AnatGrps(i).Skip))+1;
        %        Channels2Use{i} = Par.AnatGrps(i).Channels+1;
    end
else
    error('need channels number');
end

[nChannels, Channels2Use, SpatialSmoother] = DefaultArgs(varargin, { nChannels, [1:nChannels], [1 2 1]   });

BlockSize = 2^16;

InFp = fopen(FileIn, 'r');
OutFp = fopen(FileOut, 'w');
if iscell(Channels2Use)
    nShanks = length(Channels2Use);
else
    nShanks =1;
    Channels2Use = {Channels2Use};
end
while(~feof(InFp))
    Block = fread(InFp, [nChannels, BlockSize], 'short');
    CSDBlock = [];
    for s=1:nShanks
        SmoothLoss =(length(SpatialSmoother)-1);
        nCsdChannels = length(Channels2Use{s}) -2 - SmoothLoss ;
        CsdChannels{s}  = Channels2Use{s}(2+SmoothLoss/2:end-1-SmoothLoss/2);
        if length(SpatialSmoother) > 1
            CSDBlock(end+1:end+nCsdChannels,:) = -conv2(SpatialSmoother, 1 , diff(Block(Channels2Use{s},:), 2), 'valid');
        else
            CSDBlock(end+1:end+nCsdChannels,:) = -diff(Block(Channels2Use{s},:),2);
        end
    end

end
fwrite(OutFp, CSDBlock, 'short');

save([FileBase '.csd.ch'], 'CsdChannels');
fclose(InFp);
fclose(OutFp);