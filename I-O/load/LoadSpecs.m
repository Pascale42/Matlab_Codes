%function y = LoadSpecs(FileBase,State,Channels)
function [y, t, f] = LoadSpecs(FileBase,State, varargin)
[Channels, Ind] = DefaultArgs(varargin,{RepresentChan(FileBase),[]});
% if isempty(Channels)
%     Channels = ;
% end
%lsgspec=dir(['*.gspec*.' States{w}]); gspec1 = lsgspec(3).name;
load([FileBase '.gspec.' num2str(Channels(1)) '.' State],'t','f','-MAT');
if isempty(Ind) Ind = [1:length(t)]; end;
nBins = length(Ind);
y=zeros(nBins,length(f),length(Channels));
cnt=1;
for ch=Channels(:)'
    tmp =  load([FileBase '.gspec.' num2str(ch) '.' State],'y','-MAT');
    y(:,:,cnt) = log(tmp.y(Ind,:));
    cnt=cnt+1;
end

%outlier = 

%keyboard