function [ResSeg CluSeg] = LoadResSeg(FileBase,Els, Period)
%function LoadResSeg(FileBase,Els, Period)
% loads spikes from electrodes given by vector Els
% which are in segment given by Period = [Beg End} (in Par.SampleRate)
% assumes .res file is sorted in time!
nEls = length(Els);
for e=1:nEls
    ResFileName = [FileBase '.res.' num2str(Els(e))];
    CluFileName = [FileBase '.clu.' num2str(Els(e))];
    
    Nspk = 
    myRes = fscanf(fid, Nspk);
%    myRes= fread(fid, Nspk);
%    or fgetl line by line 
    

end




function [Value Position] = LocateNumPosition(FileName, Number, Order)
%locates the position of the 
NumDigits = length(num2str(Number));
match = [];
for lev=0:NumDigits
%    grepString = ['egrep -n --mmap -m 1 "^' num2str(Number) '[0-9]{' num2str(lev) '}"' FileName ];
    grepString = ['egrep -n --mmap "^' num2str(Number) '[0-9]{' num2str(lev) '}"' FileName ];
    % this effectively searches in the vicinity of Number for a spike
    % increasing the radius of vicinity 10 times every iteration of lev
    % (basically matching less and less digits of the number)
    % potenially one can make longer but more sparing search - decrease the
    % lower order digits to change the radius e.g. twice each iteration
    
    [status match] = unix(grepString);
    if ~isempty(match)
        break 
    end
end

%match = str2num(match);
switch Order
    case '<' %position of spike : t_spike<Number
        SpacePos = find(double(match)==10); % spaces correspond to end of line in grep output
        ColumnPos = find(double(match)==58); % spaces correspond to end of line in grep output
        Position = str2num(match(SpacePos(end-1)+1:ColumnPos(end)-1));
        Value = str2num(match(ColumnPos(end-1)+1:SpacePos(end)-1));
        
    case '>' % position of spike : t_spike>Number
        
        
end

return
        
        

        


