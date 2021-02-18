function h = psthplot(data,win,binwidth)

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

ymax = max(psth)*1.25;

% determine binwidth
samples = length(psth);
window = post-pre;
binwidth = window/(samples-1);

% plot
time = [pre:binwidth:post];
h = bar(time,psth,1,'w');    
axis([pre post 0 ymax])
hold on;
line([0 0],[0 ymax],'Color','k');

ylabel('spikes per second');
xlabel('time (msec)');