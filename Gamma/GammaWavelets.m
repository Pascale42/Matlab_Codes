% function GammaWavelets(FileBase, where, state,TG)
% 
% where :: regions where gamma was detected
% state :: brainstate
% TG :: LFP gamma-triggered 


function GammaWavelets(FileBase, where, state,TG)


Par = LoadPar([FileBase '.xml']);
FreqRange=[20:0.1:100];
% !!!!!!!  calculer 1 exemple pour avoir la taille de wt et avoir un freqlist

PeakPow=NaN(size(TG,2),2);
PeakFreq=NaN(size(TG,2),1);
MinFreq=NaN(size(TG,2),2);
MaxFreq=NaN(size(TG,2),2);

for g= 1:size(TG,2)
    
    % Wavelets
    [wt,f]=AnaWavelets(TG(:,g),Par.lfpSampleRate,FreqRange,'Gabor',15,0);
    disp([num2str(g) '/' num2str(size(TG,2))])
    t=[1:size(wt,1)];
    
    
    % GETTING REGION OF INTEREST
    
    Bg=median(std(abs(wt(200:1600,:)))); % get background power std (chosen freq?? 51:end)
    roi=abs(wt(600:1200,101:end)) >= 5*Bg; % get roi only for a 0.4416 sec x  35-80Hz window
          if isempty(find(roi ==1, 1));  roi=abs(wt(600:1200,101:end)) >=3*Bg;end
    ROI=zeros(size(wt,1), size(wt,2));
    ROI(600:1200,101:end)=roi; clear roi
    
%     figure('name', ['Show one wavelet result : ' num2str(g)], 'NumberTitle','off'); %
%     subplot 311; plot(TG(:,g)); axis tight
%     subplot 312; imagesc(t,f,abs(wt(:,:))'); axis xy; % coloration('dark'); waitforbuttonpress; clf;
%     subplot 313; imagesc(t, f,ROI'); axis xy
    
    
    if max(sum(ROI)) ~= 0
        
        % frequence peak
        sumROI=sum(ROI,1);
        [~,fidx]=max(sumROI); PeakFreq(g)=f(fidx);  clear fidx sumROI
        MinFreq(g,1)=f(find(sum(ROI,1) >=1, 1, 'first')); MinFreq(g,2)=find(sum(ROI,1) >=1, 1, 'first');
        MaxFreq(g,1)=f(find(sum(ROI,1) >=1, 1, 'last')); MaxFreq(g,2)=find(sum(ROI,1) >=1, 1, 'last');
        
        % puissance peak
        powfreq= zeros(size(wt));% fils a NaN matrix with only power values in ROI
        for x=1:size(wt,1)
            for y=1:size(wt,2)
                if ROI(x,y) ==1
                    powfreq(x,y)=abs(wt(x,y));
                end
            end
        end; clear x y
        
        pow=sum(powfreq(:,MinFreq(g,2):MaxFreq(g,2)),2);
        [PeakPow(g,1), PeakPow(g,2)] =max(pow);
        
        
        
%       clear Bg ROI pow powfreq wt
%         g=g+1; close
        
    end
end

save([FileBase '.' mfilename '.' where '.' state '.mat'], 'PeakPow','PeakFreq','MinFreq','MaxFreq')
