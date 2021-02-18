% [p, th0, r, logZ, k, n] = RayleighTest(spkph,clu,maxclu)
%
% Tests a set of angles for directional preference.
% spkph - phase in radians: can be
% single vector - one value returned,
% matrix - a value per column,
% single vector and clu vector of the same length of indexes from
% [1:max(clu)] - returns vector of max(clu) length - e.g. per cell
% if you want the output size to be the same for different subsets of
% spikes in the same file (e.g. in different periods) you might want to
% pass maxclu, which is max(clu) for the whole file. otherwise maxclu may
% be always different depending which cells fire in which period.
% p returns p-value, th0 mean angle, and r is radial distance
% k - circular concentraion coefficient
% output can be a single structure (if just one output) 
% see Fisher, statistical analysis of circular data, p.70.
% Anton's eddition  
%should be backwards compatible with Ken's original :))
function [p, th0, r, logZ, k, n] = RayleighTest(th,clu,maxclu)
if (nargin<3 | isempty(maxclu)) & nargin>1
    maxclu = max(clu);
end
if length(th)==0
    p=NaN;
    th0 = NaN;
    r = NaN;
    return;
end

if nargin==1

    x = exp(i*th);
    m = mean(x);

    th0 = angle(m);
    r = abs(m);
    n = size(th,1);

    z = n*r.*r;
    logZ = log(z);

    p = exp(-z).*(...
        1 + (2*z-z.*z)/(4*n) - (24*z - 132*z.*z + 76*z.^3 - 9*z.^4)/(288*n*n) ...
        );
    
    %doesn't work with a matrix - ken's function needs fix
    [mu k] = VonMisesFit(th);
    
    

elseif nargin>1
    if size(th,1)~=length(clu)
        error('th and clu have to be equal length');
    end
    if min(clu)<1
        gx = find(clu>0);
        clu=clu(gx); th = th(gx);
    end
    
%    if size(th,2)==1
%     if preserveindex
%         cnt= Accumulate(clu,1,max(clu));
%         gi = find(cnt~=0);

    x = exp(i*th);
    uClu = unique(clu);
    nclu = length(uClu);
    
    %m = Accumulate(clu,x,[max(clu) 1]);
    %cnt= Accumulate(clu,1,[max(clu) 1]);
    m = accumarray(clu,x,[maxclu 1]);
    n= accumarray(clu,1,[maxclu 1]);
   
    gi = find(n~=0);
    m(gi)=m(gi)./n(gi);
    

    th0 = angle(m);
    r = abs(m);
   
    p = NaN*ones(maxclu,1);
    k = NaN*ones(maxclu,1);
    logZ = NaN*ones(maxclu,1);
    
    z = n.*r.*r;
   
    gz =z(gi); gn=n(gi);
    logZ(gi) = log(gz);
    p(gi) = exp(-gz).*(1 + (2*gz-gz.*gz)./(4*gn) - (24*gz - 132*gz.*gz + 76*gz.^3 - 9*gz.^4)./(288*gn.*gn) );
    p = p(:);
    p(p==0)=1;
    
    %compute k
    kml =[];
    kml(gi) = BesselRatInv(r(gi));
    kml =kml(:);
    myi = find(gn<=15 & kml(gi)<2);
    k(gi(myi)) = max(kml(gi(myi)) -2./(gn(myi).*kml(gi(myi))),0);

    myi = find(gn<=15 & kml(gi)>=2);
    k(gi(myi)) = (gn(myi)-1).^3.*kml(gi(myi))./(gn(myi).^3+gn(myi));

    myi = find(gn>15);
    k(gi(myi))=kml(gi(myi));
    k=k(:);


    %     else
    %         x = exp(i*th);
    %         m = Accumulate([clu clu],x,max(clu));
%         cnt= Accumulate(clu,1,max(clu));
%     gi = cnt~=0;
%     m(gi)=m(gi)./cnt(gi);
%     
% 
%     th0 = angle(m);
%     r = abs(m);
%     n = hist(clu,[1:max(clu)])';
% 
%     z = n.*r.*r;
%     logZ = log(z);
%     z =z(gi); n=n(gi);
%     p(gi) = exp(-z).*(1 + (2*z-z.*z)./(4*n) - (24*z - 132*z.*z + 76*z.^3 - 9*z.^4)./(288*n.*n) );
%     p = p(:);
%     p(p==0)=1;
end


if nargout==1
    out = struct('p',p,'th0',th0, 'r',r, 'logZ',logZ, 'k',k, 'n',n);
    
    p = out;
end    