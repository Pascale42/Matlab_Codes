% function [y, f, phi, yerr, phierr, phloc, pow] = SpectrCohFreq(eeg, fMode, VARARGIN: FreqRange, WinLengthSec, SR, Whiten, name)
%
% fMode:  'c'  :: compute and save
%              'd'  :: display; only if saved before
%             'cd' ::compute and display, not saved
%
% DefaultArgs: FreqRange=[0.1 50]; WinLengthSec=4;
%              SR(sampling rate)=1250; Whiten=1; name(string)=[];
%
% Use 'name' to save the computation as name.SpectrCohFreq.mat
%

function [y, f, phi, yerr, phierr, phloc, pow] = SpectrCohFreq(eeg, fMode, varargin)

[FreqRange, WinLengthSec, SR, Whiten, name] = ...
    DefaultArgs(varargin,{[0.1 50], 4, 1250, 1, []});

switch fMode
    
    case 'c'
        
        WinLengthSample = 2^round(log2(WinLengthSec*SR));
        nFFT = 2*WinLengthSample;
        
        if Whiten ==1
            eeg = WhitenSignal(eeg,SR*2000,2);
        end
        
        tic
        [y, f, phi, yerr, phierr, phloc, pow]= mtptchd(eeg,[],[],nFFT,SR,WinLengthSample,[],3,'linear',[],FreqRange);
        toc
        
        save([name '.' mfilename '.mat'], 'y','f','phi','yerr','phierr','phloc','pow');
        
        
    case 'd'
        
        load([name '.' mfilename '.mat']);
        nChannels=size(y,3);
        
        
        figure('name',  [name ' - ' mfilename], 'NumberTitle','off');
        
        for Ch1=1:nChannels
            for Ch2 = Ch1:nChannels
                subplot(nChannels, nChannels, Ch1 + (Ch2-1)*nChannels);
                if(Ch1==Ch2)
                    plot(f,20*log10(abs(y(:,Ch1,Ch2))+eps),'r'); grid on; % power spectrum density
                    ylabel(['power ch ' num2str(Ch1) ', dB']); %xlabel('Frequency');
                else
                    confplot(f,y(:,Ch1,Ch2),yerr(:,Ch1,Ch2,1)); grid on; % coherence density
                    ylim([-0.1 1.1]); ylabel('coherence');
                    
                    subplot(nChannels, nChannels, Ch2 + (Ch1-1)*nChannels);
                    %                     plot(f,unwrap(phi(:,Ch1,Ch2)),'k');
                    confplot(f,(phi(:,Ch1,Ch2)/pi*180),(phierr(:,Ch1,Ch2)/pi*180),(phierr(:,Ch1,Ch2)/pi*180),'k','LineWidth',1); grid on;% phase shift
                     ylabel('phase shift'); ylim([-30 30]);
                    hold on; plot(f,angle(phloc(:,Ch1,Ch2)/pi*180),'g');
                    
                end
            end
            ForAllSubplots('set(gca, ''FontSize'', 6)')
        end
        
        
    case 'cd'
        
        WinLengthSample = 2^round(log2(WinLengthSec*SR));
        nFFT = 2*WinLengthSample;
        
        if Whiten ==1
            eeg = WhitenSignal(eeg,SR*2000,2);
        end
        
        tic
        [y, f, phi, yerr, phierr, phloc, pow]= mtptchd(eeg,[],[],nFFT,SR,WinLengthSample,[],3,'linear',[],FreqRange);
        toc
        
        nChannels=size(y,3);
        
        figure('name',  [name ' - ' mfilename], 'NumberTitle','off');
        
        for Ch1=1:nChannels
            for Ch2 = Ch1:nChannels
                subplot(nChannels, nChannels, Ch1 + (Ch2-1)*nChannels);
                if(Ch1==Ch2)
                    plot(f,20*log10(abs(y(:,Ch1,Ch2))+eps),'r'); grid on; % power spectrum density
                    ylabel(['power ch ' num2str(Ch1) ', dB']); %xlabel('Frequency');
                else
                    confplot(f,y(:,Ch1,Ch2),yerr(:,Ch1,Ch2,1)); grid on; % coherence density
                    ylim([-0.1 1.1]); ylabel('coherence');
                    
                    subplot(nChannels, nChannels, Ch2 + (Ch1-1)*nChannels);
                    %                     plot(f,unwrap(phi(:,Ch1,Ch2)),'k');
                    confplot(f,(phi(:,Ch1,Ch2)/pi*180),(phierr(:,Ch1,Ch2)/pi*180),(phierr(:,Ch1,Ch2)/pi*180),'k','LineWidth',1); grid on;% phase shift
                     ylabel('phase shift'); ylim([-30 30]);
                    hold on; plot(f,angle(phloc(:,Ch1,Ch2)/pi*180),'g');
                    
                end
            end
            ForAllSubplots('set(gca, ''FontSize'', 6)')
        end
       
        
end


%% Dependencies

function [y, f, t, phi, FStats]=mtchglong(varargin)
%function [yo, fo, to, phi, FStats]=mtchglong(x,nFFT,Fs,WinLength,nOverlap,NW,Detrend,nTapers,FreqRange);
% Multitaper Time-Frequency Cross-Spectrum (cross spectrogram)
% for long files - splits data into blockes to save memory
% function A=mtcsg(x,nFFT,Fs,WinLength,nOverlap,NW,nTapers)
% x : input time series
% nFFT = number of points of FFT to calculate (default 1024)
% Fs = sampling frequency (default 2)
% WinLength = length of moving window (default is nFFT)
% nOverlap = overlap between successive windows (default is WinLength/2)
% NW = time bandwidth parameter (e.g. 3 or 4), default 3
% nTapers = number of data tapers kept, default 2*NW -1
%
% output yo is yo(f, t)
%
% If x is a multicolumn matrix, each column will be treated as a time
% series and you'll get a matrix of cross-spectra out yo(f, t, Ch1, Ch2)
% NB they are cross-spectra not coherences. If you want coherences use
% mtcohere

% Original code by Partha Mitra - modified by Ken Harris 
% and adopted for long files and phase by Anton Sirota
% Also containing elements from specgram.m

% default arguments and that
[x,nFFT,Fs,WinLength,nOverlap,NW,Detrend,nTapers,nChannels,nSamples,nFFTChunks,winstep,select,nFreqBins,f,t] = mtparam(varargin);

% allocate memory now to avoid nasty surprises later
y=complex(zeros(nFFTChunks,nFreqBins, nChannels, nChannels)); % output array
if nargout>3
    phi=complex(zeros(nFFTChunks,nFreqBins, nChannels, nChannels));
end
nFFTChunksall= nFFTChunks;
%freemem = FreeMemory;
BlockSize = 2^8;
nBlocks = ceil(nFFTChunksall/BlockSize);
%h = waitbar(0,'Wait..');
for Block=1:nBlocks
    %   waitbar(Block/nBlocks,h);
    minChunk = 1+(Block-1)*BlockSize;
    maxChunk = min(Block*BlockSize,nFFTChunksall);
    nFFTChunks = maxChunk - minChunk+1;
    iChunks = [minChunk:maxChunk];
    Periodogram = complex(zeros(nFreqBins, nTapers, nChannels, nFFTChunks)); % intermediate FFTs
    Temp1 = complex(zeros(nFreqBins, nTapers, nFFTChunks));
    Temp2 = complex(zeros(nFreqBins, nTapers, nFFTChunks));
    Temp3 = complex(zeros(nFreqBins, nTapers, nFFTChunks));
    eJ = complex(zeros(nFreqBins, nFFTChunks));
    tmpy =complex(zeros(nFreqBins,nFFTChunks, nChannels, nChannels));
    % calculate Slepian sequences.  Tapers is a matrix of size [WinLength, nTapers]
    [Tapers, V]=dpss(WinLength,NW,nTapers, 'calc');
    % New super duper vectorized alogirthm
    % compute tapered periodogram with FFT 
    % This involves lots of wrangling with multidimensional arrays.
    
    TaperingArray = repmat(Tapers, [1 1 nChannels]);
    for j=1:nFFTChunks
        jcur = iChunks(j);
        if length(nOverlap)==1
        	Seg = [(jcur-1)*winstep+1: min((jcur-1)*winstep+WinLength,nSamples)];
            Segment = x(Seg,:);
        else
        	Seg = [nOverlap(jcur)-WinLength/2+1: min(nOverlap(jcur)+WinLength/2,nSamples)];
            Segment = x(Seg,:);
        end
                
        
        
        if (~isempty(Detrend))
            Segment = detrend(Segment, Detrend);
        end;
        SegmentsArray = permute(repmat(Segment, [1 1 nTapers]), [1 3 2]);
        TaperedSegments = TaperingArray .* SegmentsArray;
        
        fftOut = fft(TaperedSegments,nFFT);
        normfac = sqrt(2/nFFT); %to get back rms of original units
        Periodogram(:,:,:,j) = fftOut(select,:,:)*normfac; 
        % Periodogram: size  = nFreqBins, nTapers, nChannels, nFFTChunks
    end	
    if nargout>4 %ccompute fstats
        U0 = repmat(sum(Tapers(:,1:2:end)),[nFreqBins,1,nChannels,   nFFTChunks]);
        Mu = squeeze(sum(Periodogram(:,1:2:end,:,:) .* conj(U0), 2) ./  sum(abs(U0).^2, 2));
        Num = abs(Mu).^2;
        Sp = squeeze(sum(abs(Periodogram).^2,2));
        chunkFS = (nTapers-1) * Num ./ (Sp ./ sq(sum(abs(U0).^2, 2))- Num );
        %	sum(abs(Periodogram - U0.*repmat(Mu,[1,nTapers,1,1])), 2);
        FStats(iChunks, :, :)  = permute(reshape(chunkFS, [nFreqBins, nChannels, nFFTChunks]),[ 3 1, 2]);
    end
    % Now make cross-products of them to fill cross-spectrum matrix
    for Ch1 = 1:nChannels
        for Ch2 = Ch1:nChannels % don't compute cross-spectra twice
            Temp1 = reshape(Periodogram(:,:,Ch1,:), [nFreqBins,nTapers,nFFTChunks]);
            Temp2 = reshape(Periodogram(:,:,Ch2,:), [nFreqBins,nTapers,nFFTChunks]);
            Temp2 = conj(Temp2);
            Temp3 = Temp1 .* Temp2;
            eJ=sum(Temp3, 2);
            tmpy(:,:, Ch1, Ch2)= eJ/nTapers;
            
            % for off-diagonal elements copy into bottom half of matrix
            if (Ch1 ~= Ch2)
                tmpy(:,:, Ch2, Ch1) = conj(eJ) / nTapers;
            end            
            
        end
    end
    
    for Ch1 = 1:nChannels
        for Ch2 = 1:nChannels % don't compute cross-spectra twice
            
            if (Ch1 == Ch2)
                % for diagonal elements (i.e. power spectra) leave unchanged
                y(iChunks,:,Ch1, Ch2) = permute(abs(tmpy(:,:,Ch1, Ch2)),[2 1 3 4]);
            else
                % for off-diagonal elements, scale
                
                y(iChunks,:,Ch1, Ch2) = permute((abs(tmpy(:,:,Ch1, Ch2).^2) ...
                    ./ (tmpy(:,:,Ch1,Ch1) .* tmpy(:,:,Ch2,Ch2))), [2 1 3 4]);
%                if nargout>3
                    phi(iChunks,:,Ch1,Ch2) = permute(angle(tmpy(:,:,Ch1, Ch2) ...
                        ./ sqrt(tmpy(:,:,Ch1,Ch1) .* tmpy(:,:,Ch2,Ch2))), [2 1 3 4]); 
 %               end
            end
        end
    end
    
    
end
%close(h);
% we've now done the computation.  the rest of this code is stolen from
% specgram and just deals with the output stage

if nargout == 0
    % take abs, and use image to display results
    newplot;
    for Ch1=1:nChannels
        for Ch2 = Ch1:nChannels
            subplot(nChannels, nChannels, Ch1 + (Ch2-1)*nChannels);
            if Ch1==Ch2
                if length(t)==1
                    imagesc([0 1/f(2)],f,20*log10(abs(y(:,:,Ch1,Ch2))+eps)');axis xy; colormap(jet);
                else
                    imagesc(t+diff(t(1:2))/2,f,20*log10(abs(y(:,:,Ch1,Ch2))+eps)');axis xy; colormap(jet);colorbar
                end
                title(['Power specgram ' num2str(Ch1)]);
            else
                %imagesc the coherogram
                imagesc(t+diff(t(1:2))/2,f,(abs(y(:,:,Ch1,Ch2)))');axis xy; colormap(jet);colorbar
                title(['Coherogram ' num2str(Ch1)]);
                
                
                %display phaseogram
                subplot(nChannels, nChannels, Ch2 + (Ch1-1)*nChannels);
                imagesc(t+diff(t(1:2))/2,f,squeeze(phi(:,:,Ch1,Ch2))');axis xy; colormap(jet);colorbar
                title(['Phasogram ' num2str(Ch1)]);
            end
        end
    end
    xlabel('Time')
    ylabel('Frequency')
end

function varargout = DefaultArgs(Args, DefArgs)
% auxillary function to replace argument check in the beginning and def. args assigment
% sets the absent or empty values of the Args (cell array, usually varargin)
% to their default values from the cell array DefArgs. 
% Output should contain the actuall names of arguments that you use in the function

% e.g. : in function MyFunction(somearguments , varargin)
% calling [SampleRate, BinSize] = DefaultArgs(varargin, {20000, 20});
% will assign the defualt values to SampleRate and BinSize arguments if they
% are empty or absent in the varargin cell list 
% (not passed to a function or passed empty)
if isempty(Args)
    Args ={[]};
end

% if iscell(Args) & isstr(Args{1}) & length(Args)==1
%     Args = Args{1};
% end
    
if ~iscell(DefArgs)
    DefArgs = {DefArgs};
end
nDefArgs = length(DefArgs);
nInArgs = length(Args);
%out = cell(nDefArgs,1);
if (nargout~=nDefArgs)
    error('number of defaults is different from assigned');
    %keyboard
end
for n=1:nDefArgs
    
    if (n>nInArgs || isempty(Args{n}))
        varargout(n) = {DefArgs{n}}; %#ok<*CCAT1>
    else 
        varargout(n) = {Args{n}};
    end
end

function [x,nFFT,Fs,WinLength,nOverlap,NW,Detrend,nTapers,nChannels,nSamples,...
        nFFTChunks,winstep,select,nFreqBins,f,t,FreqRange] = mtparam(P)
% helper function to do argument defaults etc for mt functions

nargs = length(P);

x = P{1};
if (nargs<2 || isempty(P{2})); nFFT = 1024; else nFFT = P{2}; end;
if (nargs<3 || isempty(P{3})); Fs = 1250; else Fs = P{3}; end;
if (nargs<4 || isempty(P{4})); WinLength = nFFT; else WinLength = P{4}; end;
if (nargs<5 || isempty(P{5})); nOverlap = WinLength/2; else nOverlap = P{5}; end;
if (nargs<6 || isempty(P{6})); NW = 3; else NW = P{6}; end;
if (nargs<7 || isempty(P{7})); Detrend = ''; else Detrend = P{7}; end;
if (nargs<8 || isempty(P{8})); nTapers = 2*NW -1; else nTapers = P{8}; end;
if (nargs<9 || isempty(P{9})); FreqRange = [0 Fs/2]; else FreqRange = P{9}; end
% Now do some compuatations that are common to all spectrogram functions
if size(x,1)<size(x,2)
    x = x';
end
nChannels = size(x, 2);
nSamples = size(x,1);

if length(nOverlap)==1
    winstep = WinLength - nOverlap;
    % calculate number of FFTChunks per channel
    %remChunk = rem(nSamples-Window)
    nFFTChunks = max(1,round(((nSamples-WinLength)/winstep))); %+1  - is it ? but then get some error in the chunking in mtcsd... let's figure it later
    t = winstep*(0:(nFFTChunks-1))'/Fs;
else
    winstep = 0;
    nOverlap = nOverlap(nOverlap>WinLength/2 & nOverlap<nSamples-WinLength/2);
    nFFTChunks = length(nOverlap);
    t = nOverlap(:)/Fs; 
end 
%here is how welch.m of matlab does it:
% LminusOverlap = L-noverlap;
% xStart = 1:LminusOverlap:k*LminusOverlap;
% xEnd   = xStart+L-1;
% welch is doing k = fix((M-noverlap)./(L-noverlap)); why?
% turn this into time, using the sample frequency


% set up f and t arrays
if isreal(x)%~any(any(imag(x)))    % x purely real
	if rem(nFFT,2),    % nfft odd
		select = [1:(nFFT+1)/2];
	else
		select = [1:nFFT/2+1];
	end
	nFreqBins = length(select);
else
	select = 1:nFFT;
end
f = (select - 1)'*Fs/nFFT;
nFreqRanges = size(FreqRange,1);
%if (FreqRange(end)<Fs/2)
    if nFreqRanges==1
        select = find(f>FreqRange(1) & f<FreqRange(end));
        f = f(select);
        nFreqBins = length(select);
    else
        select=[];
        for i=1:nFreqRanges
            select=cat(1,select,find(f>FreqRange(i,1) & f<FreqRange(i,2)));
        end
        f = f(select);
        nFreqBins = length(select);
    end

function [y,f, phi, yerr, phierr, phloc,pow, normmat]=mtptchd(varargin)
%function [y, f,  phi, yerr, phierr, phloc,pow] = 
%       mtptchd(x,Res, Clu, nFFT,Fs,WinLength,nOverlap,NW,Detrend,nTapers, FreqRange,CluSubset,Segments,pval,MinSpikes);
% Multitaper coherence density for continuous and point process
% %%%%%%%%%%%%%%%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%
% x : input time series, Res/Clu - spikes 
% nFFT = number of points of FFT to calculate (default 1024)
% Fs = sampling frequency (default 2)
% WinLength = length of moving window (default is nFFT)
% nOverlap = overlap between successive windows (default is WinLength/2)
% NW = time bandwidth parameter (e.g. 3 or 4), default 3
% nTapers = number of data tapers kept, default 2*NW -1
% FreqRange =[Fi_beg Fi_end]
% ClustSubset = set of indexes to use for point process estimates (default
% all)
% 
% $$$$$$$$$$$$$$$$    Output: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% y - f x ch x ch matrix of coherence/spectra, f - freq, phi - phase shift
% ch = nChannelsEeg + length(CluSubset)
% yerr: matrix nChannelsAll x nChannelsAll x 2 - lower/upper limits. Factor
% of the spectral power for lower/upper value for spectra for now (diagonal),
% offdiagonal - coherence confidence level.
% phierr = variance of the phase (theoretical)
% Original code by Partha Mitra - modified by Ken Harris
% Also containing elements from specgram.m
% and adopted for point processes, long files and phase by Anton Sirota
% Poitnt Process spectra calculation following Jarvis&Mitra,2000 and Bijan Pesaran's code

[x, Res, Clu, nFFT,Fs,WinLength,nOverlap,NW,Detrend,nTapers,nChannels, nClu, nSamples, ...
    nFFTChunks,winstep,select,nFreqBins,f,t, FreqRange, CluSubset, Segments, pval,MinSpikes] = mtparam_pt(varargin);
winstep = WinLength - nOverlap;
nChannelsAll = nChannels + nClu;
[yo, f, normmat, pow4coh, phloc] = mtptcsd(x,Res, Clu, nFFT,Fs,WinLength,nOverlap,NW,Detrend,nTapers, FreqRange,CluSubset,Segments,pval,MinSpikes);

y = zeros(size(yo));
phi = zeros(size(yo));
pow = zeros(size(yo)); %#ok<PREALL>

%if nargout>2; phi = zeros(size(y)); end

% main loop
for Ch1 = 1:nChannelsAll
    for Ch2 = 1:nChannelsAll

        if (Ch1 == Ch2)
            % for diagonal elements (i.e. power spectra) leave unchanged
            y(:,Ch1, Ch2) = abs(yo(:,Ch1, Ch2));
           
        else
            % for off-diagonal elements, scale
            normpow = squeeze(sqrt(pow4coh(:,Ch1,Ch2,1).*pow4coh(:,Ch1,Ch2,2)));
            %if ~any(normpow==0)
            nonz = normpow~=0;
            y(nonz,Ch1, Ch2) = abs(yo(nonz,Ch1, Ch2)) ./ normpow(nonz);
            %end
            if nargout>2
%                 if ~any(normpow==0)
                    phi(nonz,Ch1,Ch2) = angle(yo(nonz,Ch1, Ch2) ./normpow(nonz));
 %                end
            end
            
        end
    end
end
pow = abs(pow4coh);

if nargout>3 % assymptotic error
    pp=1-pval/2;
    qq=1-pval;
    if length(CluSubset)==1
        numsp = length(Clu);
    else
        numsp = hist(Clu,CluSubset); %number of spikes of each cell
    end
    %coherence
    yerr = []; phierr = [];
    for Ch1=1:nChannelsAll
        %spectra
        dof = 2*normmat(Ch1,Ch1);
        if dof==0
            continue;
        end
        if Ch1>nChannels
            %finite size correction for spikes
            dof = fix(1/(1./dof + 1./(2*numsp(Ch1-nChannels))));
        end
        Qp=Schi2inv(pp,dof);
        Qq=Schi2inv(qq,dof);
        yerr(1:nFreqBins,Ch1,Ch1,1) = dof*y(:,Ch1,Ch1)/Qp;
        yerr(1:nFreqBins,Ch1,Ch1,2) = dof*y(:,Ch1,Ch1)/Qq;

        %coherence
        for Ch2=Ch1+1:nChannelsAll
            dof = 2*normmat(Ch1,Ch2);
            if Ch1<=nChannels
                dof1 = dof;
            else
                %finite size correction
                dof1=fix(2*numsp(Ch1-nChannels)*dof/(2*numsp(Ch1-nChannels)+dof));
            end
            if Ch2<=nChannels
                dof2 = dof;
            else
                dof2=fix(2*numsp(Ch2-nChannels)*dof/(2*numsp(Ch2-nChannels)+dof));
            end
            dof =min(dof1,dof2);
            yerr(1:nFreqBins,Ch1,Ch2,1) = repmat(confc(dof,pval),nFreqBins,1);
            yerr(1:nFreqBins,Ch2,Ch1,1) = yerr(1:nFreqBins,Ch1,Ch2,1);

            %now phase variance
            perfi = find((squeeze(y(:,Ch1,Ch2))-1).^2 < 10^-5 | squeeze(y(:,Ch1,Ch2))==0);
            goodi = setdiff([1:nFreqBins],perfi);
            if ~isempty(goodi)
                phierr(goodi,Ch1,Ch2) = sqrt(2./dof.*(1./(y(goodi,Ch1,Ch2).^2) - 1)); %CORRECT devide by zero
            end
            if ~isempty(perfi)
                phierr(perfi,Ch1,Ch2) = 0;
            end
            phierr(:,Ch2,Ch1) = -phierr(:,Ch1,Ch2);

        end
    end
end
% plot stuff if required

if (nargout<1)
    PlotMatrix(f,y);
end;

function confC = confc(dof,p)
%helper to compute coherence confidence values
if dof <= 2
    confC = 1;
else
    df = 1./((dof/2)-1);
    confC = sqrt(1 - p.^df);
end;
return

function x = Schi2inv(p,v)
%CHI2INV Inverse of the chi-square cumulative distribution function (cdf).
%   X = CHI2INV(P,V)  returns the inverse of the chi-square cdf with V  
%   degrees of freedom at the values in P. The chi-square cdf with V 
%   degrees of freedom, is the gamma cdf with parameters V/2 and 2.   
%
%   The size of X is the common size of P and V. A scalar input
%   functions as a constant matrix of the same size as the other input.   

%   References:
%      [1]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%      Functions", Government Printing Office, 1964, 26.4.
%      [2] E. Kreyszig, "Introductory Mathematical Statistics",
%      John Wiley, 1970, section 10.2 (page 144)


if nargin < 2, 
    error('Requires two input arguments.');
end

[errorcode, p, v] = Sdistchck(2,p,v);

if errorcode > 0
    error('Requires non-scalar arguments to match in size.');
end

% Call the gamma inverse function. 
x = Sgaminv(p,v/2,2);

% Return NaN if the degrees of freedom is not positive.
k = (v <= 0);
if any(k(:))
    x(k) = NaN;
end

function x = Sgaminv(p,a,b)
%GAMINV Inverse of the gamma cumulative distribution function (cdf).
%   X = GAMINV(P,A,B)  returns the inverse of the gamma cdf with  
%   parameters A and B, at the probabilities in P.
%
%   The size of X is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.    
%
%   GAMINV uses Newton's method to converge to the solution.

%   References:
%      [1]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%      Functions", Government Printing Office, 1964, 6.5.

%   B.A. Jones 1-12-93

if nargin<3, 
    b=1;
end

[errorcode, p, a, b] = Sdistchck(3,p,a,b);

if errorcode > 0
    error('Requires non-scalar arguments to match in size.');
end

%   Initialize X to zero.
x = zeros(size(p));

k = find(p<0 | p>1 | a <= 0 | b <= 0);
if any(k),
    tmp  = NaN;
    x(k) = tmp(ones(size(k)));
end

% The inverse cdf of 0 is 0, and the inverse cdf of 1 is 1.  
k0 = find(p == 0 & a > 0 & b > 0);
if any(k0),
    x(k0) = zeros(size(k0)); 
end

k1 = find(p == 1 & a > 0 & b > 0);
if any(k1), 
    tmp = Inf;
    x(k1) = tmp(ones(size(k1))); 
end

% Newton's Method
% Permit no more than count_limit interations.
count_limit = 100;
count = 0;

k = find(p > 0  &  p < 1 & a > 0 & b > 0);
if (~any(k(:))), return; end
pk = p(k);

% Supply a starting guess for the iteration.
%   Use a method of moments fit to the lognormal distribution. 
mn = a(k) .* b(k);
v = mn .* b(k);
temp = log(v + mn .^ 2); 
mu = 2 * log(mn) - 0.5 * temp;
sigma = -2 * log(mn) + temp;
xk = exp(Snorminv(pk,mu,sigma));

h = ones(size(pk)); 

% Break out of the iteration loop for three reasons:
%  1) the last update is very small (compared to x)
%  2) the last update is very small (compared to sqrt(eps))
%  3) There are more than 100 iterations. This should NEVER happen. 

while(any(abs(h) > sqrt(eps)*abs(xk))  &&  max(abs(h)) > sqrt(eps)    ...
                                 && count < count_limit), 
                                 
    count = count + 1;
    h = (Sgamcdf(xk,a(k),b(k)) - pk) ./ Sgampdf(xk,a(k),b(k));
    xnew = xk - h;
    % Make sure that the current guess stays greater than zero.
    % When Newton's Method suggests steps that lead to negative guesses
    % take a step 9/10ths of the way to zero:
    ksmall = find(xnew < 0);
    if any(ksmall),
        xnew(ksmall) = xk(ksmall) / 10;
        h = xk-xnew;
    end
    xk = xnew;
end


% Store the converged value in the correct place
x(k) = xk;

if count == count_limit, 
    fprintf('\nWarning: GAMINV did not converge.\n');
    str = 'The last step was:  ';
    outstr = sprintf([str,'%13.8f'],h);
    fprintf(outstr);
end

function x = Snorminv(p,mu,sigma)
%NORMINV Inverse of the normal cumulative distribution function (cdf).
%   X = NORMINV(P,MU,SIGMA) finds the inverse of the normal cdf with
%   mean, MU, and standard deviation, SIGMA.
%
%   The size of X is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.    
%
%   Default values for MU and SIGMA are 0 and 1 respectively.

%   References:
%      [1]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%      Functions", Government Printing Office, 1964, 7.1.1 and 26.2.2


if nargin < 3, 
    sigma = 1;
end

if nargin < 2;
    mu = 0;
end

[errorcode, p, mu, sigma] = Sdistchck(3,p,mu,sigma);

if errorcode > 0
    error('Requires non-scalar arguments to match in size.');
end

% Allocate space for x.
x = zeros(size(p));

% Return NaN if the arguments are outside their respective limits.
k = find(sigma <= 0 | p < 0 | p > 1 | isnan(p));
if any(k)
    tmp  = NaN;
    x(k) = tmp(ones(size(k))); 
end

% Put in the correct values when P is either 0 or 1.
k = find(p == 0);
if any(k)
    tmp  = Inf;
    x(k) = -tmp(ones(size(k)));
end

k = find(p == 1);
if any(k)
    tmp  = Inf;
    x(k) = tmp(ones(size(k))); 
end

% Compute the inverse function for the intermediate values.
k = find(p > 0  &  p < 1 & sigma > 0);
if any(k),
    x(k) = sqrt(2) * sigma(k) .* erfinv(2 * p(k) - 1) + mu(k);
end

function [x, Res, Clu, nFFT,Fs,WinLength,nOverlap,NW,Detrend,nTapers,nChannels, nClu, nSamples,...
    nFFTChunks,winstep,select,nFreqBins,f,t,FreqRange,CluSubset,Segments,pval,MinSpikes] = mtparam_pt(P)
 % helper function to do argument defaults etc for mt functions
% for point process Dims - temporal dimansion of input nSamples
nargs = length(P);
% analysis of what's the input
if (~isempty(P{1}) && length(P{1}) == length(P{2}) )
    % have only point process
    Res = P{1};
    Clu = P{2};
    x=[];
    nChannels = 0;
    nSamples = max(Res)+1;
    nInputs = 2;
elseif isempty(P{1}) && length(P{2}) == length(P{3})
   % have only point process
    Res = P{2};
    Clu = P{3};
    x=[];
    nChannels = 0;
    nSamples = max(Res)+1;
    nInputs = 3;  
elseif ~isempty(P{1}) && (length(P{2}) == length(P{3})) && ~isempty(P{2}) && length(P{2})>1
    Res = P{2};
    Clu = P{3};
    x = P{1};
    % check for column vector input
    nChannels = size(x,2);
    nSamples =  size(x,1);
    if nSamples == 1 || nSamples < nChannels
        x = x';
        nSamples = size(x,1);
        nChannels = size(x,2);
    end;
    nInputs=3;
elseif ~isempty(P{1}) && isempty(P{2}) && isempty(P{3})
    x = P{1};
    nChannels = size(x,2);
    nSamples =  size(x,1);
    Clu=[];Res=[];
    if nSamples == 1 || nSamples < nChannels
        x = x';
        nSamples = size(x,1);
        nChannels = size(x,2);
    end;
    nInputs=3;
elseif ~isempty(P{1}) && ~isempty(P{2}) && length(P{2})==1
    x = P{1};
    nChannels = size(x,2);
    nSamples =  size(x,1);
    Clu=[];Res=[];
    if nSamples == 1 || nSamples < nChannels
        x = x';
        nSamples = size(x,1);
        nChannels = size(x,2);
    end;
    nInputs=1;
    
else
    error('some of the signal inputs are wrong');
end

CluInd = unique(Clu);
if ~isempty(x) 
nT = max(size(x,1));
else
nT = max(Res);
end
P = P(nInputs+1:end); nargs = nargs - nInputs ;
if (nargs<1 || isempty(P{1})); nFFT = 1024; else nFFT = P{1}; end;
if (nargs<2 || isempty(P{2})); Fs = 2; else Fs = P{2}; end;
if (nargs<3 || isempty(P{3})); WinLength = nFFT; else WinLength = P{3}; end;
if (nargs<4 || isempty(P{4})); nOverlap = WinLength/2; else nOverlap = P{4}; end;
if (nargs<5 || isempty(P{5})); NW = 3; else NW = P{5}; end;
if (nargs<6 || isempty(P{6})); Detrend = ''; else Detrend = P{6}; end;
if (nargs<7 || isempty(P{7})); nTapers = 2*NW -1; else nTapers = P{7}; end;
if (nargs<8 || isempty(P{8})); FreqRange = [0 Fs/2]; else FreqRange = P{8}; end
if (nargs<9 || isempty(P{9})); CluSubset = CluInd; else CluSubset = P{9}; end
if (nargs<10 || isempty(P{10})); Segments = [1 nT]; else Segments = P{10}; end
if (nargs<11 || isempty(P{11})); pval=0.05; else pval = P{11}; end
%MinSpikes = 2*nTapers - theoretical recommendation from Jarvis and Mitra
if (nargs<12 || isempty(P{12})); MinSpikes=2*nTapers; else MinSpikes = P{12}; end 
% Now do some compuatations that are common to all spectrogram functions


if ~isempty(Res)
    %select spikes from Segments, if not done at the input
    [Res, ind] = SelectPeriods(Res,Segments,'d',1);
    Clu = Clu(ind);
    CluInd = unique(Clu);

    %check if there are some missing clusters in Clu and exclude them from
    %CluSubset
    CluSubset = intersect(CluSubset,CluInd);
    nClu = length(CluSubset);
    myClu = ismember(Clu, CluSubset);
    Res = Res(myClu);
    Clu = Clu(myClu);
else
    CluSubset = [];
    nClu=0;
end
%here is the time information

winstep = WinLength - nOverlap;
% calculate number of FFTChunks
nFFTChunks = round(((nSamples-WinLength)/winstep));
if nFFTChunks==0
    nFFTChunks=1;
end
% turn this into time, using the sample frequency
t = winstep*(0:(nFFTChunks-1))'/Fs;
 
% set up f and t arrays
if ~any(any(imag(x))) || isempty(x)   % x purely real
	if rem(nFFT,2),    % nfft odd
		select = [1:(nFFT+1)/2];
	else
		select = [1:nFFT/2+1];
	end
	nFreqBins = length(select);
else
	select = 1:nFFT;
end
f = (select - 1)'*Fs/nFFT;
nFreqRanges = size(FreqRange,1);
if (FreqRange(end)<Fs/2)
    if nFreqRanges==1
        select = find(f>FreqRange(1) & f<FreqRange(end));
        f = f(select);
        nFreqBins = length(select);
    else
        select=[];
        for i=1:nFreqRanges
            select=cat(1,select,find(f>FreqRange(i,1) & f<FreqRange(i,2)));
        end
        f = f(select);
        nFreqBins = length(select);
    end
end

function [y, f, normmat, pow4coh, phloc]=mtptcsd(varargin)
%function [yo, fo, normmat, pow4coh, phloc]=mtptcsd(x,Res, Clu, nFFT,Fs,WinLength,nOverlap,NW,Detrend,nTapers, FreqRange,
% CluSubset,Segments,pval,MinSpikes);
% Multitaper Hybrid Spikes-Field Cross-Spectral Density
% x : input time series, Res/Clu - spikes for the
% nFFT = number of points of FFT to calculate (default 1024)
% Fs = sampling frequency (default 2)
% WinLength = length of moving window (default is nFFT)
% nOverlap = overlap between successive windows (default is WinLength/2)
% NW = time bandwidth parameter (e.g. 3 or 4), default 3
% nTapers = number of data tapers kept, default 2*NW -1
%
% output yo is yo(f)
%
% If x is a multicolumn matrix, each column will be treated as a time
% series and you'll get a matrix of cross-spectra out yo(f, Ch1, Ch2)
% NB they are cross-spectra not coherences. If you want coherences use
% mtcohere
%
% Original code by Partha Mitra/Bijan Pesaran - modified by Ken Harris
% Also containing elements from specgram.m
% and adopted for point processes, long files and other stuff by Anton Sirota
% Poitnt Process spectra calculation following Jarvis&Mitra,2000 and Bijan Pesaran's code

% default arguments and that
[x, Res, Clu, nFFT,Fs,WinLength,nOverlap,NW,Detrend,nTapers,nChannels, nClu, nSamples,...
    nFFTChunks,winstep,select,nFreqBins,f,t,FreqRange,CluSubset,Segments,pval,MinSpikes]  = mtparam_pt(varargin);

winstep = WinLength - nOverlap;
nChannelsAll = nChannels + nClu;
clear varargin; % since that was taking up most of the memory!

% allocate memory now to avoid nasty surprises later
y=complex(zeros(nFreqBins, nChannelsAll, nChannelsAll)); % output array
pow4coh=zeros(nFreqBins, nChannelsAll, nChannelsAll,2); % output array
phloc=zeros(nFreqBins, nChannelsAll, nChannelsAll); % output array

Temp1 = complex(zeros(nFreqBins, nTapers));
Temp2 = complex(zeros(nFreqBins, nTapers));
Temp3 = complex(zeros(nFreqBins, nTapers));
eJ = complex(zeros(nFreqBins,1));

% calculate Slepian sequences.  Tapers is a matrix of size [WinLength, nTapers]
[Tapers, V]=dpss(WinLength,NW,nTapers, 'calc');

% New super duper vectorized alogirthm
% compute tapered periodogram with FFT
% This involves lots of wrangling with multidimensional arrays.

TaperingArray = repmat(Tapers, [1 1 nChannels]);
nUsedUnitsChunks = zeros(nClu,1); % count chunks that had something for each cell
nUsedEegChunks = 0;
normmat = zeros(nChannelsAll, nChannelsAll);
for j=1:nFFTChunks
    if sum(WithinRanges((j-1)*winstep+[1 WinLength],Segments))==0
        continue;
    end
    Periodogram = complex(zeros(nFreqBins, nTapers, nChannelsAll)); % intermediate FFTs
    %    ComputeCsd = 0; %a priory skip csd calculation
    ComputeCsd = zeros(nChannelsAll,1); %a priory skip csd calculation
    %continuous case
    if nChannels > 0
        Segi=(j-1)*winstep+[1:WinLength];
        %now find if this chunk is within Segments (e.g. REM)
        Segi_in = SelectPeriods(Segi,Segments,'d',1);
        if length(Segi_in)>length(Segi)/2
            Segment = zeros(WinLength,nChannels);
            indinseg = ismember(Segi,Segi_in); %find indexes of the good samples within current segment
            Segment(indinseg,:) = x(Segi_in, :);

            if (~isempty(Detrend))
                Segment = detrend(Segment, Detrend);
            end;
            SegmentsArray = permute(repmat(Segment, [1 1 nTapers]), [1 3 2]);
            TaperedSegments = TaperingArray .* SegmentsArray;

            fftOut = fft(TaperedSegments,nFFT);
            Periodogram(:,:,1:nChannels) = fftOut(select,:,:)*sqrt(2/nFFT);
            nUsedEegChunks = nUsedEegChunks + 1;
            %            ComputeCsd = 1;
            ComputeCsd(1:nChannels,1)=1;
        end
    end

    %point process data periodograms
    Segment = (j-1)*winstep+[1 WinLength];
    SegmentId = find(Res>=Segment(1) & Res<=Segment(2));
    if ~isempty(SegmentId)
        SegmentRes = Res(SegmentId) ;
        SegmentClu = Clu(SegmentId);
        if length(CluSubset)==1
            numsp = length(SegmentClu);
        else
            numsp = hist(SegmentClu,CluSubset)'; %number of spikes of each cell
        end
        if any(numsp>=MinSpikes)
            SegmentRes = SegmentRes - Segment(1) + 1;
            fftOut = PointFFT(Tapers,SegmentRes, SegmentClu, nClu, CluSubset, nFFT,Fs, MinSpikes);
            Periodogram(: , : ,  nChannels+1:nChannelsAll) = fftOut(select, : , :)*sqrt(2/nFFT);

            %    nUsedUnitsChunks = nUsedUnitsChunks +ismember(CluSubset(:),Clu(SegmentId)); %accumulate counter for each cell if it fires in the segment
            nUsedUnitsChunks = nUsedUnitsChunks + double(numsp>=MinSpikes); %accumulate counter for each cell if it fires more than MinSpikes
            ComputeCsd(nChannels+1:end) = numsp>=MinSpikes;
        end
        %        ComputeCsd = 1;
    end


    % Now make cross-products of them to fill cross-spectrum matrix
    for Ch1 = 1:nChannelsAll
        for Ch2 = Ch1:nChannelsAll % don't compute cross-spectra twice
            if ComputeCsd(Ch1)  && ComputeCsd(Ch2) % don't do csd calculation if both signals don't fall into Segments
                % here is the fix - don't compute the diagonal of y if any
                % of the components 
 
                Temp1 = squeeze(Periodogram(:,:,Ch1));
                Temp2 = squeeze(Periodogram(:,:,Ch2));
                Temp2 = conj(Temp2);
                Temp3 = Temp1 .* Temp2;
                eJ=sum(Temp3, 2);
                y(:,Ch1, Ch2)= y(:,Ch1,Ch2) + eJ;

                normmat(Ch1,Ch2) = normmat(Ch1,Ch2) + nTapers;
               % also need to accumulate the power for cross-spectra
                % normalization for unit-field case - to take into account
                % that units are not firing (enough) in many segments -
                % this can create a bias 
                pow4coh(:,Ch1, Ch2,1)= pow4coh(:,Ch1,Ch2,1) + sum(Temp1.*conj(Temp1),2);
                pow4coh(:,Ch1, Ch2,2)= pow4coh(:,Ch1,Ch2,2) + sum(Temp2.*conj(Temp2),2);
                
                %phase locking index
                phloc(:,Ch1,Ch2) = phloc(:,Ch1,Ch2) + sum(exp(i*angle(Temp1)+i*angle(Temp2)),2); %since Temp2 is conjugated
                
            end
        end
    end
end

% now fill other half of matrix with complex conjugate
for Ch1 = 1:nChannelsAll
    for Ch2 = (Ch1+1):nChannelsAll % don't compute cross-spectra twice
        normmat(Ch2,Ch1) = normmat(Ch1,Ch2);
        y(:, Ch2, Ch1) = conj(y(:,Ch1,Ch2));
        pow4coh(:, Ch2, Ch1,:) = pow4coh(:,Ch1,Ch2,:);
    end
end


normp = permute(repmat(normmat,[1,1,nFreqBins]),[3 1 2]);

%y = y ./ normp; %CORRECT devide by zero
nz = normp>0;
y(nz) = y(nz) ./ normp(nz);

nz = normp>0;
phloc(nz) = phloc(nz) ./ normp(nz);

normp = repmat(normp,[1 1 1 2]);
nz = normp>0;
pow4coh(nz) = pow4coh(nz)./normp(nz);

% we've now done the computation.  the rest of this code is stolen from
% specgram and just deals with the output stage

if nargout == 0
    % take abs, and plot results
    newplot;
    for Ch1=1:nChannelsAll
        for Ch2 = 1:nChannelsAll
            subplot(nChannelsAll, nChannelsAll, Ch1 + (Ch2-1)*nChannelsAll);
            plot(f,20*log10(abs(y(:,Ch1,Ch2))+eps));
            grid on;
            if(Ch1==Ch2)
                ylabel('psd (dB)');
            else
                ylabel('csd (dB)');
            end;
            xlabel('Frequency');
        end
    end
end

 function out = WithinRanges(x, Ranges, RangeLabel, Mode)

% WithinRanges(x, Ranges, RangeLabel, Mode)
% detects which points of the input vector lie within
% one of the ranges specified in the nx2 array ranges
% returns an array the size of x with a 1 if the corresponding
% point is in the ranges.
%
% ranges are (start1 stop1 ; start2 stop2 ; etc.)
% The ranges may be optionally labeled 1..nLabels
% in which case out is a matrix with one column per
% range label
%
% endpoint behaviour is inclusive
% (if i am right that sort() leaves equal in the order they
% were given in).
%
% if Mode is 'matrix' (default) it will give a matrix output
% with 1 if the point belongs to that range label
% if Mode is 'vector' it will give a vector, specifying the range
% of each point (gives error if any point belongs to more than 1)

% reshape x to a vector
x = x(:);

% get size info
nPoints = length(x);
nRanges = size(Ranges,1);

if nargin<3 || isempty(RangeLabel)
	RangeLabel = ones(nRanges, 1);
end
nLabels = max(RangeLabel);
if nargin<4 || isempty(Mode)
    Mode = 'matrix';
end

RangeLabel = RangeLabel(:)';
if nRanges==0
    out=zeros(size(x));
    return;
end

% check End comes after Start in each case
if any(Ranges(:,2)<Ranges(:,1))
    error('End should come after Start!');
end

% make array containing points, starts and finishes

%ToSort = [x ; Ranges(:,1) ; Ranges(:,2)];
ToSort = [Ranges(:,1) ; x ; Ranges(:,2)]; % this order means it will be inclusive
% sort it
[Sorted, Index] = sort(ToSort);

% Make delta array containing 1 for every start and -1 for every stop
% with one column for each range label

if strcmp(Mode, 'matrix')
    Delta = zeros(nPoints+2*nRanges,nLabels);
    Delta(sub2ind([nPoints+2*nRanges,nLabels],1:nRanges,RangeLabel)) = 1;
    Delta(sub2ind([nPoints+2*nRanges,nLabels],nPoints+nRanges+(1:nRanges),RangeLabel)) = -1;

    %Arrange it in order
    DeltaSorted = Delta(Index,:);

    % take cumulative sums
    Summed = cumsum(DeltaSorted);

	% and reorder back to the original order
	ReOrdered = zeros(nPoints+2*nRanges,nLabels);
	ReOrdered(Index,:) = Summed;
	
	out = ReOrdered(nRanges+1:nPoints+nRanges,:);
elseif strcmp(Mode, 'vector')
    nDelta = zeros(nPoints+2*nRanges,1);
    nDelta(1:nRanges) = 1;
    nDelta(nPoints+nRanges+(1:nRanges)) = -1;
    rDelta = zeros(nPoints+2*nRanges,1);
    rDelta(1:nRanges) = RangeLabel;
    rDelta(nPoints+nRanges+(1:nRanges)) = -RangeLabel;
    
    nDeltaSorted = nDelta(Index);
    rDeltaSorted = rDelta(Index);

    % take cumulative sums
    nSummed = cumsum(nDeltaSorted);
    rSummed = cumsum(rDeltaSorted);

	% and reorder back to the original order
	nReOrdered(Index) = nSummed;
	rReOrdered(Index) = rSummed;
    
    if any(nReOrdered(nRanges+1:nPoints+nRanges)>1)
        error('Some points belong to more than one range');
    else
        out = rReOrdered(nRanges+1:nPoints+nRanges);
    end
end
    
out = out==1;

function [y, ind] = SelectPeriods(x,Periods,varargin)
[SigType,WhereFlag, ifSquash] = DefaultArgs(varargin,{'c',1,0});
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
if ischar(Periods)
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
    
    if (nargin<3 || isempty(SigType))
        if (nTimeBins < max(Periods(:)))
            SigType = 'd';
        else
            SigType = 'c';
            if (nargout >2)
                error('too many output parameters');
%                 exit
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


if (nargin <4 || isempty(WhereFlag) )
    WhereFlag =1;
end



if (SigType == 'd')
    if ~iscell(x)
        [y, ind] = SelPerDiscr(x,Periods,WhereFlag,ifSquash);
    else
        y = cell(celli,cellj);
        ind = cell(celli,cellj);
        for ii=1:celli
            for jj=1:cellj
                [y{ii,jj}, ind{ii,jj}]= SelPerDiscr(x{ii,jj},Periods,WhereFlag,ifSquash);
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

function [y, ind]=SelPerDiscr(x,Periods, WhereFlag,SquashTime)
nPeriods = size(Periods,1);
if nargin<4 || isempty(SquashTime)
    SquashTime=0;
end
x = x(:);
%x = sort(x);
nTimeBins = max(x);
ind = [];y=[];
if SquashTime
    Shift = Periods(:,1) -[0; Periods(1:end-1,2)]-1 ; % there was evil bug here , thanks to Kenji it is fixed. damn complicated me!
    Shift = cumsum(Shift);
end 

if WhereFlag
    for p=1:nPeriods
        myi   = find(x>=Periods(p,1) & x<=Periods(p,2));
        ind = [ind; myi(:)];

        if ~SquashTime
            y = [y; x(myi) ];
        else
            y = [y; x(myi)-Shift(p) ];
        end
    end

else
    OutPeriods = [];
    if Periods(1,1)>2
        OutPeriods=[1 Periods(1,1)-1];
    end
   
    for  p=1:nPeriods-1
        OutPeriods = [OutPeriods; [Periods(p,2)+1 Periods(p+1,1)-1]];
    end
    
    if Periods(end,2)<nTimeBins
        OutPeriods = [OutPeriods; [Periods(end,2)+1 nTimeBins]];
    end

    [y, ind]=SelPerDiscr(x,OutPeriods,1,SquashTime);
    
%     ind = [ind; find(x > 1 & x < Periods(1,1))];
%     if (nPeriods>1)
%         for  p=1:nPeriods-1
%             myi=find(x > Periods(p,2) & x<Periods(p+1,1));
%             ind = [ind; myi(:)];
% 
%             if SquashTime
%                 y = [y; x(myi)];                
%             
%         end
%     end
%     ind = [ind; find(x > Periods(end,2) & x<nTimeBins)];

end
