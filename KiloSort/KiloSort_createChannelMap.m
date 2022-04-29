function KiloSort_createChannelMap(chanMap, deadChan, SamplingRate, xcoords, ycoords, kcoords)

% function KiloSort_createChannelMap(chanMap, deadChan, SamplingRate, xcoords, ycoords, kcoords)
% 
% chanMap = channels layout in row vector (base +1), even dead ones
% deadChan = Dead channels : 'none' if all channels are ok, .
% SamplingRate = sampling frequency
% xcoord, ycoord =  define the horizontal (x) and vertical (y) coordinates of the channels in row vectors. For dead channels the values won't matter. 
%                                These are in micrometer here, but the absolute scaling doesn't really matter in the algorithm. 
% kcoord =  indicates which group the channel belongs to, in a row vector. 
%                   Often, multi-shank probes or tetrodes will be organized
%                   into groups of channels that cannot possibly share spikes with the rest of the probe. 
%                   This helps the algorithm discard noisy templates shared across groups. 
%
% Saves the chanMap.mat file in the folder where the function is launched


%%  create a channel map file

connected = true(length(chanMap), 1); % i.e. which ones are not dead channels
if ~ischar(deadChan)
   connected(deadChan) = 0;
end

chanMap0ind = chanMap - 1; % indices of all channels, even dead ones

fs = SamplingRate ; % sampling frequency

save('chanMap.mat',  'chanMap','connected', 'xcoords', 'ycoords', 'kcoords', 'chanMap0ind', 'fs')
disp(['chanMap.mat file created in ' pwd])
