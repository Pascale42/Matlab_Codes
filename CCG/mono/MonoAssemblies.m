function out = MonoAssemblies(FileBase,ElLoc,varargin)

load([FileBase '.mono-' ElLoc],'-MAT');

n = size(mono.From.ElClu,1);
cnt=1;
ass = [];
conn=[];
for ii=1:n
    ass(cnt,:) = mono.From.ElClu(n,:);
    if mono.To.Type==0
        if ~isempty(ass,mono.To.ElClu(n,:))
            ass(cnt+1,:) = mono.To.ElClu(n,:);
            conn(cnt+1,cnt)
    else
        
        
    end
        
   

