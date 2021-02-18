function MakeEeg(FileBase,varargin)

Par =LoadPar([FileBase '.xml']);
rs = Par.SampleRate/1250;
cmd = ['downsample ' FileBase '.dat ' FileBase '.eeg ' num2str(Par.nChannels) ' ' num2str(rs) ' ' num2str(Par.Offset)];
fp = fopen('doit','w');
fprintf(fp,'%s\n',cmd);
system('chmod +x doit');
fclose(fp);