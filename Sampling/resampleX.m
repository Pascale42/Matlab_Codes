  
function Y = resampleX(X,alpha,r)
%  
%  By: Ron Abileah, Vista Research Inc
%
%  Original version:  November 10, 2005
%  Version 1.1        December 1, 2005 
%                         - Corrected comments and mentioned similarity to
%                         MATLAB function "resample."
%                         - Output Y is row (column) if input X is row
%                         (column) 
%                         - Sets alpha to default value (1) if alpha is not provided
%  Version 1.2        December 15, 2005
%                         - Corrected indexing bug discovered by Eike
%                         Rietsch.  This improved accuracy.
%
%  Resamples X(n).  Y(n) = X(alpha*n), where alpha is a resample interval.
%  For example, if X is data sampled at 1000 samples per second and you 
%  would like to transform it to the equivalent of 1100 samples per second
%  use alpha= 1000/1100 (.9091); for 800 sample per second 
%  use alpha = 1000/800 (1.25).
%
%  ResampleX is similar to the MATLAB "resample" function (in the Signal
%  Processcing Toolbox).  There are two differences:
%
%  (1) The MATLAB resample does some fancy schmancy interpolation of the 
%       original time series;  resampleX works on the Fourier trasnform of
%       the time series.  The main benefit of FT processing is speed.   
%  (2) MATLAB resample cputime depends on the value of alpha.  It runs faster with 
%       simple rational numbers.  The cputime of resampleX is indpendent of
%       alpha.   
%
%  The main reason for using resampleX instead of resample is speed.  In
%  test cases resampleX was generaly 5-20 times faster.    Use resample if
%  your alpha values are simple rational numbers or numerical accuracy is more
%  important.  Use resampleX for very general values of alpha or where some
%  accuracy can be traded for speed.
%
%  ResampleX uses the fact that resampling a time series X by a factor alpha is 
%  equivalent to resampling the frequency samples of its transform by m/alpha, 
%  where m is a frequency index. 
%
%  The calling sequence is one of the following
%
%             Y = resampleX(X,alpha)
%             Y = resampleX(X,alpha,r)
%             Y = resampleX
%
%  X can be real or complex. Output Y is the same length as X.  
%  If alpha > 1 some Y's will be extrapolations beyond the end of X.
%  Extrapolated values are not realiable, so throw them away.  Keep only 
%  the first N/alpha values, where N is the original length of the data.
%  Whne X is real, Y's may have small imaginary values due to
%  approximations.  
%
%  Optional parameter r is an integer frequency interpolation factor.  The
%  function will work with r = 1, but r = 2,4, or more produces more
%  accurate results.  Use r=64 if you want very accurate results and are
%  not too concerned about computing time.  The default value is r = 8.
%
%  Calling the function with no arguments produces a test signal and a plot
%  of the original and resampled signal with alpha =0.95, r = 8.  
%  

%  -- Set r to its default value if not specified as input argument
if ~exist('r')
    r = 8; end

%  -- Produce test signal if there are no input arguments
if nargin == 0 
    x1=0:pi/64:2*pi;
    X=exp(i*x1)+0.5*exp(i*3*x1)+0.25*exp(i*x1.^2);     
    alpha = 0.95;                                                     
end

% The default value of alpha is 1
if nargin == 1
    alpha =1;
end

% -- zero pad X for frequency interpolation 
Y=X;
N0=length(Y);  N=r*N0; N2 = N/2;
if r >1
    Y(N)=0;
end

%  -- Fourier transform the padded time series
Y = fft(Y); Y(N2+1)=0;  %  FFT and set -Nyquist to zero
n = 1:length(X);

%  -- Resample the Fourier transform
m = round((0:(N2-1))/alpha) - (0:(N2-1)); 
f = floor(alpha*N2); m(f:end)=NaN;
m = [ m  0 -fliplr(m(2:end)) ] + (1:N) ; 
m(find(isnan(m)))=N2+1; 
Y=ifft(Y(m));                               
Y=Y(1:N0)/alpha;
 

if nargin>0
    return
end

%  -- Display results 
figure 
n=1:N0; 
plot( n ,real(X),'-k' ,n ,imag(X),'--k',...
    n ,real(Y),'-r' ,n ,imag(Y),'--r')
xlabel('Sample number'); ylabel('Amplitude')
axis([1 N0 -2 2])
legend ('Original real', 'Original imaginary ', 'Resampled real','Resampled imaginary') 

end
