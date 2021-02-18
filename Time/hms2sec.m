% function sec=hms2sec(h,m,s)
% 
% Converts times input in hours-minutes-seconds (h m s) format to the equivalent measure in seconds.

function sec=hms2sec(h,m,s)

sec=s+(m*60)+(h*3600);
disp([num2str(sec) ' seconds'])
