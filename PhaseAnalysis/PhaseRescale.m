%function SpkPh = PhaseRescale(Ph,Res,Scales,Shifts,Win)
% Ph is from -pi:pi ...it rescales the phase Ph by variable parameter and
% get's the rescaled phase values of spikes for each scale value
% output: SpkPh = as many columns of Res length as scales

function SpkPh = PhaseRescale(Ph,Res,Scales,Shifts,Win)
nScales = length(Scales);
nShifts = length(Shifts);
%SpkPh = zeros(length(Res), nScales);

if Win>0
len = length(Ph);
nlen = len - rem(len,Win);
Ph = reshape(Ph(1:nlen),Win,[]);
end
for ii=1:nScales
    for jj=1:nShifts
        %rescale Ph vector
        ScPh = mod(Scales(ii)*unwrap(Ph)-pi,2*pi)-pi;
        ScPh = ScPh +Shifts(jj);
        if Win>0
            ScPh = reshape(ScPh,[],1);
            ScPh(end+1:end+rem(len,125))=zeros(rem(len,nlen),1);
        end
        SpkPh(:,ii,jj) = ScPh(Res);
    end
end
