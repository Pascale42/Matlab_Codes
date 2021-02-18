%function [y, ind] = SelectPeriods(x,Periods,SigType,WhereFlag, ifSquash)
%
% ex :: [eegSta orind]= SelectPeriods(eeg(:),STA,'c',1);
%
% selects from the signal (time series  - 'd', continuous - 'c')
% values in/out the Periods ranges (should be the same sampl.rate
% WhereFlag defines in(1) or out(0) (default - 1 = in)
% ifSquash is 1 if you apply it for discrete signal and want to squash the
% gaps to match the sample indexes of continuous signal
%   output: y - new signal, ind - indexes of the old signal
function [y, ind] = SelectPeriods(x,Periods,varargin)
[SigType,WhereFlag, ifSquash] = DefaultArgs(varargin,{'c',1,0});

if isstr(Periods)
    Periods = load(Periods);
end
if isempty(Periods)
    y = x;
    ind = [1:length(x)];
    return
end


Periods(find(Periods(:)==0))=1;

nPeriods  = size(Periods,1);

if ~iscell(x)
    nChannels = min(size(x));
    if size(x,1)==nChannels
        x=x';
    end
    nTimeBins = max(size(x));
    
    if (nargin<3 | isempty(SigType))
        if (nTimeBins < max(Periods(:)))
            SigType = 'd';
        else
            SigType = 'c';
            if (nargout >2)
                error('too many output parameters');
                exit;
            end
        end
    end
    if (SigType =='d')
        Channels = size(x,1);
    end
else
    celli=size(x,1);
    cellj=size(x,2);
end


if (nargin <4 | isempty(WhereFlag) )
    WhereFlag =1;
end



if (SigType == 'd')
    if ~iscell(x)
        [y ind] = SelPerDiscr(x,Periods,WhereFlag,ifSquash);
    else
        y = cell(celli,cellj);
        ind = cell(celli,cellj);
        for ii=1:celli
            for jj=1:cellj
                [y{ii,jj} ind{ii,jj}]= SelPerDiscr(x{ii,jj},Periods,WhereFlag,ifSquash);
            end
        end
    end
    
    
end


if (SigType == 'c')
    y=[];
    ind=[];
    
    
    if WhereFlag
        for p=1:nPeriods   
            y = [y; x(Periods(p,1):Periods(p,2),:)];
            ind = [ind; [Periods(p,1):Periods(p,2)]'];
        end                
    else
        if Periods(1,1)>2
            y = [y; x(1:Periods(1,1)-1,:)];
            ind = [ind; [1:Periods(1,1)-1]'];
        end
        for p=1:nPeriods-1
            y = [y; x(Periods(p,2)+1:Periods(p+1,1)-1,:)];       
            ind = [ind; [Periods(p,2)+1:Periods(p+1,1)-1]'];
        end
        
        y = [y; x(Periods(nPeriods,2)+1:end,:)];
        ind = [ind; [Periods(nPeriods,2)+1:size(x,1)]'];
        
    end
    
    %     if (nargout >2)
    %         error('too many output parameters');
    %         ind =[];
    %     end
end

