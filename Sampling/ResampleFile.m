%ResampleFile(inname,outname,numchannel,resampldown, resamplup)
% resampling program --- requires signal processing toolbox
% Hajime Hirase ... no warranty ..(sorry) (1999)
% function deciall(inname,outname,numchannel,resampl)
% this function passes anti aliasing low pass filter first
% to the data and then resamples with 1/resampl sampling rate)

function ResampleFile(inname,outname,numchannel,resampldown, resamplup)

if nargin <4,
error('function ResampleFile(inname,outname,numchannel,resampldown, resamplup)');
return
end
if nargin<5 | isempty(resamplup) %length(resampl)==1
	resamplup = 1;
end 


% open input file and output file
datafile = fopen(inname,'r');
outfile = fopen(outname,'w');
%
buffersize = 2^22-mod(2^22, resampldown);
overlaporder   = 8;%lcm(8,resamplup);
overlaporder2  = overlaporder/2*resamplup;
overlaporder21 = overlaporder2+1;
obufsize = overlaporder * resampldown;
obufsize11 = obufsize - 1;

% the first buffer

[obuf,count] = fread(datafile,[numchannel,obufsize],'int16');
obuf = fliplr(obuf);
frewind(datafile);
[datasegment,count] = fread(datafile,[numchannel,buffersize],'int16');  

%append overlap chunk to avoid edge effect
datasegment2 = [obuf,datasegment]';
resampled =resample(datasegment2,resamplup,resampldown);

%the idea here is to save only middle portion unaffected by the edge
%effects
count2 = fwrite(outfile,resampled(overlaporder2*2+1:size(resampled,1)-overlaporder2,:)','int16');
obuf = datasegment2(size(datasegment2,1)-obufsize11:size(datasegment2,1),:);
% do the rest

while ~feof(datafile),
  [datasegment,count] = fread(datafile,[numchannel,buffersize],'int16');  
  datasegment2 = [obuf;datasegment'];
  resampled = resample(datasegment2,resamplup,resampldown);
  count2 = fwrite(outfile,resampled(overlaporder21:size(resampled,1)-overlaporder2,:)','int16');
  obuf = datasegment2(size(datasegment2,1)-obufsize11:size(datasegment2,1),:);
end  

% add the last unprocessed portion 
resampled = resample(obuf,resamplup,resampldown);
count2 = fwrite(outfile,resampled(overlaporder21:end,:)','int16');
fclose(outfile);

