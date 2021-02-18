%function FileCSDAllCh(FileBase, nChannels, Channels2Use, SpatialSmoother, KeepNotUsed)
% does a CSD on a file, creates a file .csd  
%smoothes with kernel [1 2 1] by default
function FileCSDAllCh(FileBase, varargin)
FileIn = [FileBase '.eeg'];
FileOut = [FileBase '.csd'];
if FileExists([FileBase '.xml'])
    Par = LoadPar([FileBase '.xml']);
    nChannels = Par.nChannels;
    nShanks = length(Par.AnatGrps);
    for i=1:nShanks
%        if length(Par.AnatGrps(i).Channels)<3; continue; end
%        Channels2Use{i} = Par.AnatGrps(i).Channels(find(~Par.AnatGrps(i).Skip))+1;
        Channels2Use{i} = Par.AnatGrps(i).Channels+1;
    end
else
    error('need channels number');
end

[nChannels, Channels2Use, SpatialSmoother, KeepNotUsed] = ...
    DefaultArgs(varargin, { nChannels, Channels2Use, [1 2 1] , 1  });

BlockSize = 2^14;

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
%    Block = detrend(Block,'constant');

    if KeepNotUsed
        CSDBlock = Block;
    else
        CSDBlock = zeros(size(Block));
    end
    Block = ButFilter(Block',2,[0.1 250]/Par.lfpSampleRate*2,'bandpass')';
    for s=1:nShanks
        SmoothLoss =(length(SpatialSmoother)-1);
        nCsdChannels = length(Channels2Use{s}) -2 - SmoothLoss ; %number of channels where csd is non-zero (meaningfull)
        nchInShank = length(Channels2Use{s});
        csdchind = [2+SmoothLoss/2:nchInShank-1-SmoothLoss/2];
        CsdChannels{s}  = Channels2Use{s}(csdchind);
        tmpcsd=zeros(nchInShank,size(Block,2));
        if nchInShank<3; continue; end
        if length(SpatialSmoother) > 1
            tmpcsd(csdchind,:) = -conv2(SpatialSmoother, 1 , diff(Block(Channels2Use{s},:), 2), 'valid');
        else
            tmpcsd(csdchind,:) = -diff(Block(Channels2Use{s},:),2);
        end
        if KeepNotUsed
             CSDBlock(Channels2Use{s}(csdchind),:) =  tmpcsd(csdchind,:);
        else
             CSDBlock(Channels2Use{s},:) =  tmpcsd;
        end
    end
    fwrite(OutFp, int16(CSDBlock), 'short');
end
save([FileBase '.csd.ch'], 'CsdChannels');
fclose(InFp);
fclose(OutFp);

