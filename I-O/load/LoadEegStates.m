%
% function LoadEegStates(FileBase, State, fMode, Channel)
%
% !!!! fMode =  'get' or 'make'
% Channel =  channel ref number or 'tout'


function LoadEegStates(FileBase, State, fMode, Channel);



switch fMode
    
    case 'make'
        
        Par = LoadPar([FileBase '.xml']);
        STA = load([FileBase '.sts.' State]);
        
        if isnumeric(Channel);
            eeg  = LoadBinary([FileBase '.eeg'], Channel, Par.nChannels);
            [eeg orind]= SelectPeriods(eeg(:),STA,'c',1);
            Ref= Channel;
            save([FileBase '.eeg.' State '.' Ref '.mat'], 'eeg', 'orind','Par', 'Ref');
        end
        
        if isstr(Channel);
            Channel = length(Par.AnatGrps(1).Channels);
            eeg  = LoadBinary([FileBase '.eeg'], Channel, Par.nChannels);
            [eeg orind]= SelectPeriods(eeg(:),STA,'c',1);
            save([FileBase '.eeg.' State '.' Channel '.mat'], 'eeg', 'orind','Par', 'Ref');
        end
        
        
        
    case 'get'
        
        load([FileBase '.eeg.' State '.' Ref '.mat'])
end
