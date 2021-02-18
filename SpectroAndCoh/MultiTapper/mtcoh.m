function mtcoh(x, 
nt = size(x,2);
nch = size(x,2);

if nargin < 4 sampling = 1; end 
nt = nt./sampling;
if nargin < 3 tapers = [nt, 5, 9]; end 
if length(tapers) == 2
    n = tapers(1);
    w = tapers(2);
    p = n*w;
    k = floor(2*p-1);
    tapers = [n,p,k];
    disp(['Using ' num2str(k) ' tapers.']);
end
if length(tapers) == 3  
    tapers(1) = tapers(1).*sampling;  
    tapers = dpsschk(tapers); 
end
if nargin < 5 fk = [0,sampling./2]; end
if length(fk) == 1
    fk = [0,fk];
end
if nargin < 6 pad = 2; end
if nargin < 7 flag = 11; end 
if nargin < 8 pval = 0.05;  end

N = length(tapers(:,1));
if N ~= nt*sampling 
    error('Error: Tapers and time series are not the same length'); 
end

K = length(tapers(1,:));
nf = max(256,pad*2.^(nextpow2(N+1)));
nfk = floor(fk./sampling.*nf);

% Determine outputs
f = linspace(fk(1),fk(2),diff(nfk));
errorchk = 0;
if nargout > 4 errorchk = 1; end

coh = zeros(1, diff(nfk));
S_X = zeros(1, diff(nfk));
S_Y = zeros(1, diff(nfk));

coh_err = zeros(2, diff(nfk));
SX_err = zeros(2, diff(nfk));
SY_err = zeros(2, diff(nfk));

Xk = zeros(nch*K, diff(nfk));
Yk = zeros(nch*K, diff(nfk));
for ch=1:nch
    tmp1 = detrend(X(ch,:))';
    xk = fft(tapers(:,1:K).*tmp1(:,ones(1,K)),nf)';
    Xk((ch-1)*K+1:ch*K,:) = xk(:,nfk(1)+1:nfk(2));
    tmp2 = detrend(Y(ch,:))';
    yk = fft(tapers(:,1:K).*tmp2(:,ones(1,K)),nf)';
    Yk((ch-1)*K+1:ch*K,:) = yk(:,nfk(1)+1:nfk(2));
end
S_X = mean(abs(Xk).^2);
S_Y = mean(abs(Yk).^2);
coh = mean(Xk.*conj(Yk))./sqrt(S_X.*S_Y);
if errorchk			%  Estimate error bars using Jacknife
    for ik = 1:nch*K
        indices = setdiff([1:nch*K],[ik]);
        Xj = Xk(indices,:);
        Yj = Yk(indices,:);
        jcoh(ik,:)=sqrt(2*nch*K-2)*atan(Xj.*conj(Yj) ...
            ./sqrt(mean(abs(Xj).^2).*mean(abs(Yj).^2)));
        jXlsp(ik,:) = log(mean(abs(Xj).^2,1));
        jYlsp(ik,:) = log(mean(abs(Yj).^2,1));
    end
end