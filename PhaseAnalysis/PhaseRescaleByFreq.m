%function SpkPh = PhaseRescaleByFreq(Ph,Res,Freqs,Shifts)
% Ph is from -pi:pi ...it rescales the phase Ph by variable parameter and
% get's the rescaled phase values of spikes for each scale value
% output: SpkPh = as many columns of Res length as scales

function SpkPh = PhaseRescale(Ph,Res,Freqs,Shifts)
win =125;
nFreqs = length(Freqs);
nShifts = length(Shifts);
%SpkPh = zeros(length(Res), nScales);

%make segments out of continuyous phase
len = length(Ph);
nlen = len - rem(len,win);

Ph = reshape(Ph(1:nlen),win,[]);

uPh =unwrap(Ph);
%fit lsqr lines to unrwapped phases to get mean freqs. on each
%segment
for kk=1:size(Ph,2)
    pol(kk,:)=polyfit([1:win]',uPh(:,kk),1);
end
mFr = pol(:,1);
%norm the unwrapped phase by the freq.
nPh = uPh./repmat(mFr(:)',win,1);
for ii=1:nFreqs
    for jj=1:nShifts
        %rescale Ph vector
        ScPh = mod(Freqs(ii)*nPh+Shifts(jj)-pi,2*pi)-pi;
        ScPh = reshape(ScPh,[],1);
        ScPh(end+1:end+rem(len,win))=zeros(rem(len,nlen),1);
        SpkPh(:,ii,jj) = ScPh(Res);
    end
end
