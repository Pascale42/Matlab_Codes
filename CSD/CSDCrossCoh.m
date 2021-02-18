function CSDCrossCoh(filename, Channels, FreqBand, csdStep, Periods, IfIn)


ifFilter = 0;
%fn = 'sp804.013';
%Freqband =[8 18];
csd =[];
csd = [csd MakeCSD(filename,Channels,csdStep,ifFilter)];
nch=size(csd,2);
coh = zeros(nch,nch);
phi = zeros(nch,nch);
for i=1:nch
    for j=i+1:nch
        [c f p] = mtchglong(csd(:,[i j]),2^8,1250,2^10,2^9,[],'linear',[],FreqBand);
       %keyboard
        [coh(i,j) ind]= max(c(:,1,2));
        phi(i,j) = p(ind,1,2);
    end
end

%rem =load('rem');
%rcsd = SelectedPeriods(csd,        