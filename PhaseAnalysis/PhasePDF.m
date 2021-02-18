function [PDF, f, PhBins] = PhasePDF(varargin)
%Phase Prob. Density Function
% code is stolen from TFE function
%   PDF = PhasePDF(X,Y,NFFT,Fs,WINDOW) estimates the transfer function of the 
%   system with input X and output Y using Welch's averaged periodogram 
%   method.  X and Y are divided into overlapping sections, each of which 
%   is detrended, then windowed by the WINDOW parameter, then zero-padded 
%   to length NFFT.  The magnitude squared of the length NFFT DFTs of the 
%   sections of X are averaged to form Pxx, the Power Spectral Density of X.
%   The products of the length NFFT DFTs of the sections of X and Y are 
%   averaged to form Pxy, the Cross Spectral Density of X and Y.  Txy
%   is the quotient of Pxy and Pxx; it has length NFFT/2+1 for NFFT even, 
%   (NFFT+1)/2 for NFFT odd, or NFFT if X or Y is complex. If you specify 
%   a scalar for WINDOW, a Hanning window of that length is used.  Fs is 
%   the sampling frequency which does not effect the transfer function 
%   estimate but is used for scaling of plots.
%
%   [Txy,F] = TFE(X,Y,NFFT,Fs,WINDOW,NOVERLAP) returns a vector of freq-
%   uencies the same size as Txy at which the transfer function is 
%   estimated, and overlaps the sections of X and Y by NOVERLAP samples.
%
%   PhasePDF(X,Y,...,DFLAG), where DFLAG can be 'linear', 'mean' or 'none', 
%   specifies a detrending mode for the prewindowed sections of X and Y.
%   DFLAG can take the place of any parameter in the parameter list
%   (besides X and Y) as long as it is last, e.g. TFE(X,Y,'mean');
%   


PhMin = -2*pi; PhMax = 2*pi; nPhBins =100;
error(nargchk(2,7,nargin))
x = varargin{1};
y = varargin{2};
dflag = 'none';
nfft= varargin{3};
Fs = varargin{4}; 
window = hanning(varargin{5});
if nargin<6 | isempty(varargin{6})
    noverlap = varargin{5}/2;
end
% [msg,nfft,Fs,window,noverlap,p,dflag]=psdchk(varargin(3:end),x,y);
% error(msg)
    
% compute PSD and CSD
window = window(:);
n = length(x);		% Number of data points
nwind = length(window); % length of window
if n < nwind    % zero-pad x , y if length is less than the window length
    x(nwind)=0;
    y(nwind)=0;  
    n=nwind;
end
x = x(:);		% Make sure x is a column vector
y = y(:);		% Make sure y is a column vector
k = fix((n-noverlap)/(nwind-noverlap))	% Number of windows
					% (k = fix(n/nwind) for noverlap=0)
index = 1:nwind;

% Pxx = zeros(nfft,1); Pxx2 = zeros(nfft,1);
% Pxy = zeros(nfft,1); Pxy2 = zeros(nfft,1);
PDFtot=zeros(nfft,k);
for i=1:k
    if strcmp(dflag,'none')
        xw = window.*x(index);
        yw = window.*y(index);
    elseif strcmp(dflag,'linear')
        xw = window.*detrend(x(index));
        yw = window.*detrend(y(index));
    else
        xw = window.*detrend(x(index),0);
        yw = window.*detrend(y(index),0);
    end
    index = index + (nwind - noverlap);
    Xx = fft(xw,nfft);
    Yy = fft(yw,nfft);
%     PhDif =  phase(Yy) - phase(Xx);
    PhDif =  angle(Yy) - angle(Xx);
   %[pdfHist PhBins]= histp(PhDif',PhMin, PhMax, nPhBins);
    PDFtot(:,i) =  PhDif;%pdfHist';
     %keyboard   
%     Xx2 = abs(Xx).^2;
%     Xy2 = Yy.*conj(Xx);
%     Pxx = Pxx + Xx2;
%     Pxx2 = Pxx2 + abs(Xx2).^2;
%     Pxy = Pxy + Xy2;
%     Pxy2 = Pxy2 + Xy2.*conj(Xy2);
end

% posPhind = find(PDFtot>PhMax);
% negPhind = find(PDFtot<PhMin);
% PDFtot(posPhind)=mod(posPhind,PhMax);
% PDFtot(negPhind)=mod(negPhind,PhMin);



[PDF PhBins]= histp(PDFtot',PhMin, PhMax, nPhBins);
PDF= PDF';
% Select first half
if ~any(any(imag([x y])~=0)),   % if x and y are not complex
    if rem(nfft,2),    % nfft odd
        select = [1:(nfft+1)/2];
    else
        select = [1:nfft/2+1];   % include DC AND Nyquist
    end
    PDF = PDF(select,:);
else
    select = 1:nfft;% compute PSD and CSD

end
PDF = PDF/k;
f= (select - 1)'*Fs/nfft;

% phl = -170*pi/180;phr = 170*pi/180;phrange=10*pi/180;
% 
% fi = find(f>5 &f<11);
% lind = find(PDFtot(fi,select)>(phl-phrange) & PDFtot(fi,select)<(phl+phrange));
% rind = find(PDFtot(fi,select)>(phr-phrange) & PDFtot(fi,select)<(phr+phrange));
% 
% figure
% t=[lind; rind];
% g = [ones(length(lind),1); 2*ones(length(rind),1)];
% CCG(t,g,1,100,Fs/noverlap,unique(g),'hz');
% 
% keyboard
% set up output parameters
PhBins = PhBins*180/pi;



if nargout==0
  
    PDF=conv2(PDF,ones(3,6)/18);
    imagesc(PhBins,f,PDF);
    ylabel('Frequency'), xlabel('Phase');
    ylim([0 100]);
    colorbar
end



function [num, bins] = histp(x, xmin, xmax, nbins)
%HISTP	Histogram estimate of 1-dimensional probability distribution.
%
%	Description
%
%	HISTP(X, XMIN, XMAX, NBINS) takes a column vector X  of data values
%	and generates a normalized histogram plot of the  distribution. The
%	histogram has NBINS bins lying in the range XMIN to XMAX.
%
%	H = HISTP(...) returns a vector of patch handles.
%
%	See also
%	DEMGAUSS
%
%	Copyright (c) Ian T Nabney (1996-9)
%   Anton, modified for many columns
ndata = size(x,1);
bins = linspace(xmin, xmax, nbins);
binwidth = (xmax - xmin)/nbins;
num = hist(x, bins);
num = num/(ndata*binwidth);
return