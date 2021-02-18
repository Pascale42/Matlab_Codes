% function CleanEegFrames(FileBase,   )
%
% Must create a .sts.Clean file before
% FR = frame rate
% Output = textfile of beginning/end frames to cut out on AviDemux software
%


function CleanEegFrames(FileBase, FR)

if ~FileExists([FileBase '.sts.Clean'])
    disp('Run CheckEegStates and Create FileBase ".sts.Clean" first!!!')
    return
end


% load the Par and the Clean segments
Par=LoadPar([FileBase '.xml']);
clean=load([FileBase '.sts.Clean']);
clean=clean/Par.lfpSampleRate;

% Get the frames to cut out
frames=clean*FR; % video @ 25 fps
frames=round(frames);
framestocut=NaN(size(frames, 1)+1, size(frames,2));

for n=1:size(frames, 1)
    framestocut(n, 2)=frames(n, 1);
    framestocut(n+1, 1)=frames(n, 2);
end

% Write in a text file
fid=fopen([FileBase '.clean.txt'], 'w'); 
fprintf(fid,'%s','Frames To Cut Out (beginning/end frames)'); fprintf(fid,'\n');
for n = 1:size(framestocut,1) %
fprintf(fid,'%g\t',framestocut(n,:)); fprintf(fid,'\n');
 end 
fclose(fid);
