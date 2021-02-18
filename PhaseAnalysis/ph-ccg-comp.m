
ccg = CatStruct(out,{'ccg','t','phhist','phbin'},4);
t=sq(ccg.t(1,:,1,1));

figure(121);clf
for i=1:29
subplotfit(i,29)
bar(t,sq(ccg.ccg(:,1,2,i)));axis tight
xlim([-150 150]);
end

phbin=sq(ccg.phbin(1,:,1,1));
figure(122);clf
for i=1:29
subplotfit(i,29)
phh= sq(ccg.phhist(1,:,1,i));
bar([phbin phbin+360],[phh phh]);axis tight
%xlim([-150 150]);
end


%ForAllSubplots('axis tight');