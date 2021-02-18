function h = rasterplot1(data,win);

SPIKEHEIGHT = 1;

pre = win(1);
post = win(2);

j=length(data);
w=length(data{1,1});

for(i = 2:j)
    w = max(w,length(data{1,j}));
end

sr = j;
sc = w;

y = [1:sr]';
Y = repmat(y,1,sc);

spikes = zeros(sr,sc);
for(i =1:j)
    spikes(i,1:length(data{1,i}))=data{1,i};
end

N = sr*sc;
SpikeVec(1:N) = spikes(1:N);
TrialVec(1:N) = Y(1:N);
blanks = find(isnan(SpikeVec));
SpikeVec(blanks) = [];
TrialVec(blanks) = [];

s = SpikeVec;
y0 = TrialVec;
y1 = TrialVec-SPIKEHEIGHT;
h = line([s;s],[y0;y1],'Color','k');
axis([pre post 0 sr]);
% axis off;
hold on;
line([0 0],[0 sr],'Color','k');

xlabel('time (msec)');
ylabel('trial');