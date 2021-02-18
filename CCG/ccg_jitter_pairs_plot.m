%
% function ccg_jitter_pairs_plot(FileBase, pairs)


function ccg_jitter_pairs_plot(FileBase, pairs)

for n=133:size(pairs,1)

CCG_jitter_one(FileBase,pairs(n,[1 2]),pairs(n,[3 4]),'d');

 print(gcf, '-depsc', [FileBase '.CCGJitterOne_' num2str(pairs(n,1)) '.'  num2str(pairs(n,2)) '-' num2str(pairs(n,3)) '.' num2str(pairs(n,4)) '.eps'] );
end
