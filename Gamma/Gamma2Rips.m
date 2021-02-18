% function Gamma2Rips(FileBase, State)

function Gamma2Rips(FileBase, State)

%
% switch fMode
%
%     case 'compute'


Par = LoadPar([FileBase '.xml']);
STA = load([FileBase '.sts.' State]);

% loads all channels

load([FileBase '.DetectGammaBursts.mPFC.' State '.mat'], 'GBtime');
gam = GBtime; % for a given state


if ~FileExists([FileBase '.spw.mat'])
    fprintf('no spws');
    return;
end
load([FileBase '.spw.mat']);
rips = Rips.t;

figure
scatter(1:length(rips), rips)
hold on
scatter(1:length(gam), gam, 'xm')