function imagesclogy(times,freqs,data,clim, xticks, yticks, varargin)

  if size(data,1) ~= length(freqs)
      fprintf('imagesclogy(): data matrix must have %d rows!\n',length(freqs));
      return
  end
  if size(data,2) ~= length(times)
      fprintf('imagesclogy(): data matrix must have %d columns!\n',length(times));
      return
  end
  if min(freqs)<= 0
      fprintf('imagesclogy(): frequencies must be > 0!\n');
      return
  end
  
  steplog = log(freqs(2))-log(freqs(1)); % same for all points
  realborders = [exp(log(freqs(1))-steplog/2) exp(log(freqs(end))+steplog/2)];
  newfreqs = linspace(realborders(1), realborders(2),length(freqs));
  
  % regressing 3 times
  border = mean(newfreqs(2:end)-newfreqs(1:end-1))/2; % automatically added to the borders in imagesc
  newfreqs = linspace(realborders(1)+border, realborders(2)-border,length(freqs));
  border = mean(newfreqs(2:end)-newfreqs(1:end-1))/2; % automatically added to the borders in imagesc
  newfreqs = linspace(realborders(1)+border, realborders(2)-border,length(freqs));
  border = mean(newfreqs(2:end)-newfreqs(1:end-1))/2; % automatically added to the borders in imagesc
  newfreqs = linspace(realborders(1)+border, realborders(2)-border,length(freqs));
  
  if nargin == 4 & ~isempty(clim)
      imagesc(times,newfreqs,data,clim);
  else 
      imagesc(times,newfreqs,data);
  end;
  set(gca, 'yscale', 'log');
  
  % puting ticks
  % ------------
  if nargin >= 5, set(gca, 'xtick', xticks); end;
  if nargin >= 6
      divs = yticks;
  else 
      divs = linspace(log(freqs(1)), log(freqs(end)), 10);
      divs = ceil(exp(divs)); divs = unique(divs); % ceil is critical here, round might misalign out-of border label with within border ticks
  end;
  set(gca, 'ytickmode', 'manual');
  set(gca, 'ytick', divs);
  
  % additional properties
  % ---------------------
  set(gca, 'yminortick', 'off', 'xaxislocation', 'bottom', 'box', 'off', 'ticklength', [0.03 0], 'tickdir','out', 'color', 'none');
  if ~isempty(varargin)
      set(gca, varargin{:});
  end;