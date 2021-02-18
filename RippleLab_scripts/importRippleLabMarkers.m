% function importRippleLabMarkers(structure, method, State)
%
% Loads gamma markers from RippleLab
%
% INPUTS
% structure :: string of the structure name where gamma was detected
% method   ::  'beginning' for importing only beginning of events; 'all' for beginnings and ends
% State       :: 'SWS'...
%
% OUTPUT
% FileBase.DetectGammaBursts.structure.mat, GBtime

function GBtime=importRippleLabMarkers(structure, varargin)



[method, State] = DefaultArgs(varargin,{'beginning', 'SWS'});

FileBase=gfb;

if FileExists([FileBase '.' State '_SLL' structure '.mat'])
    OutArgs=load([FileBase '.' State '_SLL' structure '.mat']);
else
    disp([[FileBase '.' State '_SLL' structure '.mat'] ' does not exist!!'])
    return
end


switch method
    case 'beginning'
        
        GBtime = OutArgs.(genvarname([structure])).st_HFOInfo.m_EvtLims;
        GBtime=GBtime(:,1);
        

            
            if FileExists([FileBase '.' State '_VIS' structure '.mat'])
                OutArgs2=load([FileBase '.' State '_VIS' structure '.mat']);
                GBtime2 = OutArgs2.(genvarname([structure])).st_HFOInfo.m_EvtLims;
                GBtime2=GBtime2(:,1);
                
                GBtime=sort([GBtime; GBtime2]);
            else
                disp([([FileBase '.' State '_VIS' structure '.mat']) ' does not exist: No VIS detection added'])
            end

        


 case 'all'
        
        GBtime = OutArgs.(genvarname([structure])).st_HFOInfo.m_EvtLims;
        

            
            if FileExists([FileBase '.' State '_VIS' structure '.mat'])
                OutArgs2=load([FileBase '.' State '_VIS' structure '.mat']);
                GBtime2 = OutArgs2.(genvarname([structure])).st_HFOInfo.m_EvtLims;
                
                GBtime=sort([GBtime; GBtime2]);
            else
                disp([([FileBase '.' State '_VIS' structure '.mat']) ' does not exist: No VIS detection added'])
            end

        
end

% save([FileBase '.DetectGammaBursts.' structure '.mat'], 'GBtime');