function psthmean = psthsetup(data,win,binwidth);

pre = win(1);
post = win(2);


time = pre-binwidth:binwidth:post-binwidth;

j=length(data);
H=zeros(j,length(time));

for(i = 1:j)
    spikes = data{1,i};
    H(i,:) = histc(spikes',time)';
end

H = H*(1000/binwidth);
psth = mean(H);

psth(1) = [];
psth(length(psth)) = [];

psthmean=psth;