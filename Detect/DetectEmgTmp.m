%function DetectEmg(FileBase,Regime,Overwrite)
% Regime: 1 - do the detection, 0 - just compute info. and store
function DetectEmg(FileBase,varargin)

[Regime,Overwrite]  = DefaultArgs(varargin,{1,0});

FileName = [FileBase '.DetectEmg.mat'];

SampleRate =1250;
Par = LoadPar([FileBase '.par']);
EmgPar = LoadEmgPar(FileBase);

if Regime==0 | ~FileExists(FileName) %do the computation

	emg = bload([FileBase '.emg.dat'],[EmgPar.nChannels inf]);
	emg= emg';
	
	
	
	%WinInSeconds = 0.1; % seconds
	%optWi		ndow = 2^(round(log2(SampleRate*WinInSeconds))+1);
	%[spec, f, t ] = mtcsglong(WhitenSignal(emg), optWindow*2, SampleRate, optWindow,optWindow/2,[],[],[],[1 400]);
	
	%for i=1:EmgPar.nChannels
	%    subplot(EmgPar.nChannels,1,i);
	%    imagesc(t,f,log(spec(:,:,i))');
	%end
	%pow = sq(mean(log(spec),2));
	%dpow = sq(diff(pow,1,1));
	
	%powPeak = 
	
	
	
	
	%PlotTraces(emg(1:1250*20,:));
	%SlowEmg = ButFilter(emg,2,[3 40]/(SampleRate/2), 'bandpass');
	FastEmg = ButFilter(emg,2,[65 300]/(SampleRate/2), 'bandpass');
	
	% Sigma = 0.1;
	SmoothLen = SampleRate*20/1000;
	% SmoothWin = [-SmoothLen:SmoothLen]/SmoothLen;
	% Smoother = exp(-SmoothWin.^2/Sigma^2/2);
	% Smoother = Smoother./sum(Smoother);
	
	%FastEmg = SmoothMatrix(FastEmg,SmoothWin);
	
	%do resampling
	SamplCoef =10; % RESAMPLING COEFFICIEN~!!!!!!!!!!!!!!!!!!!!!!!!
	FastEmg = resample(FastEmg,1,SamplCoef);
	SampleRate = SampleRate/SamplCoef;
	
	FastEmg = abs(FastEmg);	
	FastEmg = Filter0(ones(SmoothLen,1)/SmoothLen,FastEmg);
	
	%rf = resample(FastEmg,1,10);
	%rs = resample(abs(SlowEmg),1,10);
	
	%now get peaks on each channel
	
	EmgRes = [];EmgClu = [];
	MaxInt = 50*SampleRate/1000;
	BurstThr = round(500*SampleRate/1000); 
	BurstDuration = round(30*SampleRate/1000);
	for i=1:EmgPar.nChannels
		th = 3*std(sq(FastEmg(:,i)));
		myPks = LocalMinima(-FastEmg(:,i),MaxInt, -th);
		EmgRes = [EmgRes; myPks];
		EmgClu = [EmgClu; ones(length(myPks),1)*i];
	end
	[dummy sortind ] = sortrows([EmgRes(:) EmgClu(:)],1);
	EmgRes = EmgRes(sortind);
	EmgClu = EmgClu(sortind);
	
	[EmgB.Burst, EmgB.BurstLen, EmgB.SpkPos, EmgB.OutOf] = SplitIntoBursts(EmgRes,BurstThr);
	EmgBurst = EmgRes(EmgB.Burst);
	
	EmgSegs = GetSegs(FastEmg, EmgBurst - BurstThr, 2*BurstThr);
	
	EmgBefore = sq(sum(EmgSegs(1:BurstThr-BurstDuration,:,:),1));
	EmgAfter = sq(sum(EmgSegs(BurstThr+BurstDuration:end,:,:),1));
	
	EmgPow = sq(sum(EmgSegs(BurstThr-BurstDuration:BurstThr+BurstDuration,:,:),1));
	EmgRes = EmgRes*SamplCoef;
	EmgBurst = EmgBurst*SamplCoef;
	
	save(FileName,'EmgRes','EmgClu','EmgBurst','EmgBefore','EmgAfter','EmgPow','EmgB');
	clear emg
end

if Regime==1 % do the manual control 
	load(FileName);
	Fet.Len = EmgB.BurstLen;
	nch = length(EmgPar.Channels);
	uniqueLabels = unique(EmgPar.Labels);
	nlab = length(uniqueLabels);
	
	for i=1:nlab
		labind = find(strcmp(EmgPar.Labels,uniqueLabels{i}));
		Fet.(['Pow_' uniqueLabels{i}]) = mean(EmgPow(:,labind),2);
	end
	Fet.Before_All =  mean(EmgBefore,2);
	Fet.After_All =  mean(EmgAfter,2);
	
	for i=1:nlab
		labind = find(strcmp(EmgPar.Labels,uniqueLabels{i}));
		Fet.(['Before_' uniqueLabels{i}]) = mean(EmgBefore(:,labind),2);
	end
	for i=1:nlab
		labind = find(strcmp(EmgPar.Labels,uniqueLabels{i}));
		Fet.(['After_' uniqueLabels{i}]) = mean(EmgAfter(:,labind),2);
	end
	%Fet = [BurstLen, EmgPow, EmgBefore, EmgAfter];
	%Clu = KlustaKwik(Fet,5,10);
	%keyboard
	[Clu,dummy,GoodInd] = xgobi(Fet);
	Clu = Clu(:); GoodInd = GoodInd(:);
	%keyboard
	
	% save into events file
	MakeEvtFile(EmgBurst(GoodInd),[FileBase '.emg.evt'],Clu,SampleRate);
	msave([FileBase '.emgb'],[EmgBurst(GoodInd)*16 Clu]);
end








%-------------------------------------------
%Dat = readmulti([FileBase '.eeg'],Par.nChannels, Channels);
% 
% [icasig, A, W] = fastica(Dat,'displayMode', 'on');
% nDims = size(A,2);
% BadICAs = input('Enter ICAs that correspond to artifacts?');
% GoodICAs = setdiff([1:nDims],BadICAs);
% RecDat = A(:,GoodICAs)*icasig(GoodICAs,:);
% 
% bsave([FileBase '.emg'],[Dat' ; RecDat]);

%WhiteRecDat = WhitenSignal(Dat);



