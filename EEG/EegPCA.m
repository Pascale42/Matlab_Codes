%data = readsingle('l23-01-018.eeg.0')';
function EegPCA(data,FRange)

fftl= 2^9;

n_segment = floor(length(data)'/fftl);

sr = 1250;

n_freq = fftl/2 + 1;

f_axis = [0:sr/fftl:sr/2];
axislimit = [FRange(1):FRange(2)];
axislimit2 = axislimit;

Spectra = zeros(n_segment,n_freq);

for jj=1:n_segment
  data_begin = (jj-1)*fftl+1;
  data_end = jj*fftl;
  datapiece = data(data_begin:data_end);
  datapiecep = abs(fft(hanning(fftl).*datapiece)).^2;
  Spectra(jj,:) = datapiecep(1:n_freq)';
end

Clean_Spectra = Spectra;

Selected_Spectra = Clean_Spectra(:,axislimit);
Selected_Spectra2 = Clean_Spectra(:,axislimit2);
Spectra_Mean = mean(Selected_Spectra);
Spectra_Std  = std(Selected_Spectra);

C_Spectra = Selected_Spectra - repmat(Spectra_Mean,n_segment,1);
N_Spectra = C_Spectra ./ repmat(Spectra_Std,n_segment,1);

[pc,score,latent] = princomp(N_Spectra);

Projected = N_Spectra * pc;

figure;
clf

subplot(8,2,1);
plot(f_axis(axislimit),pc(:,1));
title('Principal Components');
subplot(8,2,3);
plot(f_axis(axislimit),pc(:,2));
subplot(8,2,5);
plot(f_axis(axislimit),pc(:,3));
subplot(8,2,7);
plot(f_axis(axislimit),pc(:,4));
xlabel('frequency');


subplot(4,2,2);
imagesc([1,n_segment]*fftl/sr,...
	[f_axis(axislimit2(1)),f_axis(axislimit2(end))],...
	log(Selected_Spectra2)');
axis xy;
xlabel('time (s)');
ylabel('freq (Hz)');
title('Specgram');

subplot(4,2,4);
imagesc([1,n_segment]*fftl/sr,...
	[f_axis(axislimit(1)),f_axis(axislimit(end))],...
	N_Spectra');
axis xy;
xlabel('time (s)');
ylabel('freq (Hz)');
title('Norm. Specgram');

subplot(2,2,3);
plot(Projected(:,1),Projected(:,2),'.','MarkerSize',1);
xlabel('1st PC');ylabel('2nd PC');

subplot(4,2,6);

CProjected1 = conv (ones(5,1)/5,Projected(:,1));
SProjected1 = CProjected1(3:end-2);

plot([1:n_segment]*fftl/sr,SProjected1,'r-');
hold on
plot([1:n_segment]*fftl/sr,Projected(:,1),':');


xlabel('time (s)');
ylabel('1st PC');

subplot(4,2,8);
plot([1:n_segment]*fftl/sr,Projected(:,2));
xlabel('time (s)');
ylabel('2nd PC');
letterlayout
