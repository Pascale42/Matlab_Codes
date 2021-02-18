% function prteps(fname,location,varargin)
% [plus] = DefaultArgs(varargin,{[]});
% fname - name of the file
% location = 1  saves it in current directory
% location = '+' save it in '/plus/' directory in the current location 
% or in any location '/...../'
%

function prteps(fname,location,varargin)
 [plus] = DefaultArgs(varargin,{[]});

if location == 1
  location= [pwd '/']; 
end

if strcmp(location, '+') == 1
  location=[pwd plus];
end

print(gcf, '-depsc2','-painters',['-r600'],[location fname]);
disp(['Saves ' fname ])

% 
% function prteps(fname,varargin)
% [figh,resol,location] = DefaultArgs(varargin,{gcf,600,1});
% if location==1
%     upath = '/ant3/antsiro';
% elseif isstr(location)
%     upath = location;
% else
%     upath = userpath;
%     sl = strfind(upath,'/');
%     upath = upath(1:sl(end)-1);
% end
% print(figh, '-depsc2','-painters',['-r' num2str(resol)],[upath '/figures/theta/' fname]);
% %print(figh, '-depsc2','-zbuffer',['-r' num2str(resol)],[upath
% %'/figures/theta/' fname]);
