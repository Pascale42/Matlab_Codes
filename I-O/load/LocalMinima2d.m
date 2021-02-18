function out = LocalMinima2d(x)

cent = round(nT/2);
[nT nCh] = size(x);

bigch = x(cent,:)>mean(x,1)+std(x,0,1);
bigch = find(bigch);
x = x(:,bigch); %now only analyze bigch
%hoping that it will capture the right minima
[minamp minti] = min(x);

[mintis newchind] = sort(minti);

%find the first minima (channelwise)
[min1 firsti] = min(minti);
minti = minti - minti(firsti);

%now find which channel takes minimal value
[min2 minchi] = min(mint);

out.Amp = zeros(nCh,1)*NaN;
out.Lat = zeros(nCh,1)*NaN;

%out.Chans = bigch;%(newchind);
out.Amp(bigch) = minamp;
out.Lat(bigch) = minti;




