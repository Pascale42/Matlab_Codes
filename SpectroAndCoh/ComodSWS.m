%function [Co, fo] = ComodSWS(filename, Ch, nFFT, Fs, WinLength, nOverlap,FreqRange, NW, Detrend, nTapers)
function out = ComodSWS(filename, varargin)
shuffle =1;
[Ch, nFFT, Fs,WinLength, nOverlap, FreqRange, NW, Detrend, nTapers] = DefaultArgs(varargin, ...
{[], 2^10,  1250, 2^10, 0, [0 300], 3, [], 5});

x= cxhipeeg(filename, Ch);
rem = load([filename '.rem']);
rem = round(rem/16);
x = SelectPeriods(x,rem,'c',0);
for shift = [0 2];
if shuffle
    x = [x(:,1) [x(shift*WinLength+1:end,2); x(1:shift*WinLength,2)]];
end

Clip = 0;
% if Clip is 1, correlations below 0 will be replaced by 0.

nChannels = size(x,2);
nSamples = size(x,1);

% calculate spectrograms
[y, f, t]=mtchglong(x,nFFT,Fs,WinLength,nOverlap,NW,Detrend,nTapers,FreqRange);
nFreqBins = size(y,2);
nTimeBins = size(y,1);

spex =[];
% FreqBins = find(f >= FreqRange(1) & f <= FreqRange(2));
% fo = f(FreqBins);
% main loop
for i=1:nChannels
			ynorm = y(:, :, i, i);
			ynorm = ynorm ./ repmat(mean(ynorm,1),  nTimeBins,1);
            spex(:,:,i) = log(ynorm);
end

%spex = abs(spex);

% calculate correlation coefficients
DataMat = reshape(spex, [nTimeBins, nFreqBins*nChannels]);
				
CorrMat = corrcoef(DataMat);
if (Clip) CorrMat = clip(CorrMat, 0, 1); end;

% produce output array and plot(if required)
C = zeros(nFreqBins,nFreqBins);

for i=1:nChannels
	for j=1:nChannels
		
% 		C(:,:) = CorrMat((i-1)*nFreqBins + (1:nFreqBins), (j-1)*nFreqBins + (1:nFreqBins));
		C = CorrMat((i-1)*nFreqBins + (1:nFreqBins), (j-1)*nFreqBins + (1:nFreqBins));
        if (i==1 & j==2)
            out{1}(:,:,1) = C;
            if (shift ==0)
                out{1}(:,:,2) = C;
            else
                out{1}(:,:,2) = out{1}(:,:,2) -C;
            end
        end
  		if (nargout<1)
 			subplot(nChannels, nChannels, j + (i-1) * nChannels);
			imagesc(f, f, C(:,:));
            set(gca,'ydir','norm');
        end;
		
		drawnow;
	end;
end;
end
out{2}= f;