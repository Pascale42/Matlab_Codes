DbName = 'mice';
%change here to get other type
Restrict = {'s','', '', '', 0 , 1};
db = RestrictDb(DbName,  Restrict);

o2s = GoThroughDb('Osc2Spikes2P',db, 0, 1, 1, 5, 30, 5, 30, [1 2 3], [1:6],1);
% o2s = GoThroughDb('Osc2Spikes2',db, 0, 1, 1, 10, 50, 5, 50, [1 7 9], []);
% 
% m=[];
% for i=1:8
% for ref = 1:3
% %    ImageMatrix(o2s{i}{3}, o2s{i}{2}, squeeze(o2s{i}{6}(: ,: ,ref ,: ,:)));
% %    pause
%     m(:,:,ref,1:6,1:6,i) = squeeze(o2s{i}{6}(: ,: ,ref ,1:6 ,1:6));
% end
% end
% 
dccg = c2m(o2s,1);
rccg = c2m(o2s,4);
xccg = c2m(o2s,5);

mdccg = mean(dccg,6);
mrccg = mean(rccg,5);
mxccg = mean(xccg,5);

tr = o2s{1}{2};
tx = o2s{1}{3};

for i=1:3
%     figure
%     PlotMatrix(tr,squeeze(mrccg(:,i,:,:)));
    
    figure
    PlotMatrix(tx,squeeze(mxccg(:,i,:,:)));
    
%     figure
%     ImageMatrix(tr,tx,squeeze(mdccg(:,:,i,:,:)));
end
