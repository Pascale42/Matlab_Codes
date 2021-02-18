

% function CleanEegTS(FileBase)
%
% Must create a .sts.Clean file before


function CleanEegTS(FileBase)

if ~FileExists([FileBase '.sts.Clean'])
    disp('Run CheckEegStates and Create FileBase ".sts.Clean" first!!!')
    return
end


% load the Par and the Clean segments
Par=LoadPar([FileBase '.xml']);
clean=load([FileBase '.sts.Clean']);
clean=clean/Par.lfpSampleRate;

% Get the episodes to cut out

tocut=NaN(size(clean, 1)+1, size(clean,2));

for n=1:size(clean, 1)
    tocut(n, 2)=clean(n, 1);
    tocut(n+1, 1)=clean(n, 2);
end

% Write in a text file
fid=fopen([FileBase '.tocut.txt'], 'w'); 
fprintf(fid,'%s','Frames To Cut Out (beginning/end frames)'); fprintf(fid,'\n');
for n = 1:size(tocut,1) %
fprintf(fid,'%g\t',tocut(n,:)); fprintf(fid,'\n');
 end 
fclose(fid);
