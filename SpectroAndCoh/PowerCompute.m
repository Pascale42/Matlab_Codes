%function PowerCompute(filename,Chans, FreqRange, Overwirte, nChannels, Window)
% computes the wavelet transformations and gives a power
function [pow, wav] = PowerCompute(filename,varargin)
[Chans, FreqRange, Overwrite, nChannels, Window] = DefaultArgs(varargin,{{'c'} , [6 18; 20 50], 0, 0,50});

global memcont;
if FileExists([filename '.par'])
    par = LoadPar([filename '.par']);
    if ~nChannels
        nChannels = par.nChannels;
    end
else
    warning('par file doesnot exist');
end
fn =[filename '.osc.mat'];
if FileExists(fn) & ~Overwrite
    return;
end

if iscell(Chans)
    Chans = GetChannels(filename, Chans,0);
end
if isempty(Chans)
    fprintf('no channels to analyze');
    return
end
% if nargout>0 & FileExists(fn)
%     savedChans = load(fn,'-MAT','Chans');
%     savedFreqRange = load(fn,'-MAT','FreqRange');
%     %    if isequal(Chans, savedChans) ...
%     %           & isequal(FreqRange, savedFreqRange)
%     if ~Overwrite
%         %pow = load(fn,'-MAT','pow');
%         %wav = load(fn,'-MAT','wav');
%         nTotSignals = sum(cellfun('length',savedFreqRange)/2);
%         pow = bload([filename '.oscpow'],[nTotSignals inf], 0, 'int16');
%         wav = bload([filename '.oscwav'],[nTotSignals inf], 0, 'int16');
%         pow = pow'; wav = wav';
%         return
%         
%     end
% end

par = LoadPar([filename '.par']);
sr = round(1e6/par.SampleTime/8);
dname = pwd;
age = strfind(pwd,'p');
if ~isempty(age)
    age = str2num(dname(age+1));
else
    age =0;
end
if age<1 and age>7
    fprintf('not a baby');
    return
end
eeg = readmulti([filename '.eeg'],nChannels,Chans);

nSignals = length(Chans);
% for many channels selected FreqRange is cell array
% do some checks of the input 
if nSignals>1 
    if ~iscell(FreqRange)
        if ndims(FreqRange)==2
            warning('assume FreqRange is the same for all channels selected');
            for k=1:nSignals
                newFreqRange{k} =  FreqRange;
            end
            FreqRange = newFreqRange;
        elseif ndims(FreqRange)==3
            if  size(FreqRange,3)==nSignals
                for k=1:nSignals
                    newFreqRange{k} =  FreqRange(:,:,k);
                end 
            else
                warning('number of freq ranges does not match signals number');
                return
            end
        end
    else
        if length(FreqRange) ~=nSignals
            warning('number of freq ranges does not match signals number');
            return
        end
    end
else
    newFreqRange{1} = FreqRange;
    FreqRange = newFreqRange;
end
nTotSignals = sum(cellfun('length',FreqRange)/2);
%eeg = resample(eeg,1,5); 
%rem = load([filename '.rem']);
%rem = round(rem/16);
%[eeg, ind] = SelectPeriods(eeg,rem,'c',0);


%now do the actual wavelete filtering 
nTime = size(eeg,1);
pow=[];wav=[];
win = sr*Window;
nWins = ceil(nTime/win);
h=waitbar(0,'Wait..');
cnt=1;
maxpow = [];
wav = zeros(nTime, nTotSignals);
pow = zeros(nTime, nTotSignals);
for i=1:nWins
    segwin = (1+(i-1)*win):min(i*win,nTime);
    for j=1:nSignals
        thisFreqRange = FreqRange{j};
        nFreqRange = size(thisFreqRange,1);
        memcont(end+1) =FreeMemory;
        [x, w] = WaveletPowerRec(eeg(segwin,j),thisFreqRange,sr, 6,1);
        memcont(end+1) =FreeMemory;
        %        fprintf('\b\b\b%d %d',i,j);
        wav(segwin,[cnt:(cnt+nFreqRange-1)]) =  w;
        pow(segwin,[cnt:(cnt+nFreqRange-1)]) =  x;
        cnt= cnt+nFreqRange;
    end
    if mod(i,50)==0
        keyboard
    end
    waitbar(i/nWins,h);
end
close(h);

if FileExists(fn) & Overwrite
    system(['mv ' filename '.oscpow ' filename '.oscpow.old']);
    system(['mv ' filename '.oscwav ' filename '.oscwav.old']);
    fidpow = fopen([filename '.oscpow'],'a+');
    fidwav = fopen([filename '.oscwav'],'a+');
else
    fidpow = fopen([filename '.oscpow'],'a+');
    fidwav = fopen([filename '.oscwav'],'a+');
end
fwrite(fidpow, pow', 'uint32');
fwrite(fidwav, wav', 'int16');

fclose(fidpow); fclose(fidwav);
if ~FileExists(fn) | Overwrite
    system(['mv ' fn ' ' fn '.old']);
    save(fn,'FreqRange', 'Chans');
end
if nargout>0
    warning('not ready yet');
end


