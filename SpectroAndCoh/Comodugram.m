% [Co, f, Pval, CoRlo, CoRup] = Comodugram(x, nFFT, SampleRate, FreqRange, WinLength, NW, Detrend, Clip)
%
% Takes an input sequence x and does a multi-pane
% plot showing correlated changes in power in frequency
% bands 
%
% addition Pascale 12/06/2015
% OUTPUT: 
% Co :: correlation coeeficient matrix
% P   :: matrix of p-values for testing the hypothesis of no correlation.  Each p-value is the probability
%   of getting a correlation as large as the observed value by random chance, when the true correlation is zero.  If P(i,j) is small, say
%   less than 0.05, then the correlation R(i,j) is significant.
% CoRlo, CoRup ::  matrices of the same size as R, containing lower and upper bounds for a 95%
%   confidence interval for each coefficient.
%
%
% nFFT and SampleRate are just like for specgram  function
%
% FreqRange = [fLow fHigh] allows you to view only a certain 
% frequency range.  Specify it in Hz.If 2 rows - freq.range is different 
% for each signal 
%
% NW is an argument for the multitaper method, as is Detrend
%
% if Clip is 1, correlations below 0 will be replaced by 0.
%
% optional output argument nTimeBins gives number of points in the regression
% - so you can do significance testing.
%
% Example: [Co nTimeBins] = Comodugram(eeg([1,5],:)',1024,1250,[1 100]);
%
%         Where: eeg is the filename loaded with bload, 1 and 5 are the channels
%                you've selected to compare (the : here is where you would put
%                a time restriction, such as a REM period). 1024 is nFFT, 1250
%                is the sampling rate, and [1 100] is the frequency range you
%                wish to compare.
%
%          NOTE: you MUST transpose the first input statement (x)

function [Co, f, Pval, CoRlo, CoRup] = Comodugram(x, varargin)

[nFFT, SampleRate, FreqRange, WinLength, NW, Detrend, Clip] = DefaultArgs(varargin, ...
    {2^11, 1250, [0 100], 2^10, 3, [], 0});




nChannels = size(x,2);
nSamples = size(x,1);

% calculate spectrograms

[spex,f,t] = mtcsglong(x, nFFT, SampleRate,WinLength,[],NW,Detrend,[],FreqRange);
spex = log(2*(NW*2-1)*spex)  - log(repmat(mean(spex), [size(spex,1),1, 1]));
spex = abs(spex);

% find frequency bins to consider

nTimeBins = size(spex,1);
nFreqBins = size(spex,2);

% calculate correlation coefficients
DataMat = reshape(spex, [nTimeBins, nFreqBins*nChannels]);
[CorrMat ,pval,rlow,rup]= corrcoef(DataMat);% CorrMat = corrcoef(DataMat); modif Pascale le 12/02/2015
      
if (Clip); CorrMat = clip(CorrMat, 0, 1); end;

% produce output array and plot(if required)
C = zeros(nFreqBins,nFreqBins); %#ok<*PREALL>
P = zeros(nFreqBins,nFreqBins);
Rlow = zeros(nFreqBins,nFreqBins);
Rup = zeros(nFreqBins,nFreqBins);


for i=1:nChannels
    for j=1:nChannels
        
        C = CorrMat((i-1)*nFreqBins + (1:nFreqBins), (j-1)*nFreqBins + (1:nFreqBins));
        P = pval((i-1)*nFreqBins + (1:nFreqBins), (j-1)*nFreqBins + (1:nFreqBins));
        Rlow = rlow((i-1)*nFreqBins + (1:nFreqBins), (j-1)*nFreqBins + (1:nFreqBins));
        Rup = rup((i-1)*nFreqBins + (1:nFreqBins), (j-1)*nFreqBins + (1:nFreqBins));
        
        if (nargout<1)
            subplot(nChannels, nChannels, j + (i-1) * nChannels);
            imagesc(f, f, C(:,:)); set(gca,'ydir','norm'); drawnow;
            
            
            
        else
            Co(:,:,i,j) = C(:,:);
            Pval(:,:,i,j) = P(:,:);
            CoRlo(:,:,i,j) = Rlow(:,:);
            CoRup(:,:,i,j) = Rup(:,:);
        end
   
    end
end


% 
% cm=flipud(jet);
%      colormap(cm);