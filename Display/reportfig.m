%function reportfig(FigHandle, FileName, Preview, Comment, Resolution,SaveFig)
% This is a quick and simple function to generate an HTML
% page containing the figure and comments on it. Figures entered
% under the same filename (or on the same day) will be added to the same
% page. The html file and directory with figures are stored in 
% homedir/mrep. Both jpeg and matlab .fig are stored. 
% FigHandle  (default = gcf) 
% FileName  - (default = current date) just name, no path 
% Preview (default = 0) if 1 - launches browser to preview the page
% Comment - to put below, default - Figure # 
function reportfig(varargin)
 
TodayDate = date;
DefFileName = ['mrep.' TodayDate];	

[FigHandle, FileName, Preview,Comment,Resolution,SaveFig] = ...
    DefaultArgs(varargin,{gcf,DefFileName,0,[],100,1});
UserPath = getenv('HOME');
LongFileName = [UserPath '/mrep/' FileName '.html'];
DirName = [UserPath '/mrep/' FileName];
if exist([UserPath '/mrep'])~=7
	mkdir(UserPath,'mrep');
end
if FileExists(LongFileName)
	fprintf('adding to existing report\n');
	DirCont = dir(DirName);
	DirCont = DirCont(3:end);
	maxindex=1;
	for i=1:length(DirCont)
		[d1,name,ext] = fileparts(DirCont(i).name);
        digind = strfind(name,'-')+1;
        digind = digind(end);
		maxindex = max(maxindex, str2num(name(digind:end)));
	end
	FigIndex = maxindex+1;
%	keyboard
else
	mkdir([UserPath '/mrep'],FileName);
        FigIndex = 1;
end
LongFigName = [DirName '/' FileName '-' num2str(FigIndex)];
FigName = [FileName '-' num2str(FigIndex)];
fPointer=fopen(LongFileName,'at+');
if fPointer<3
	error('Can''t open html file');
end		
if FigIndex == 1
	fprintf(fPointer,'<html>\n<head>\n<title>Matlab report %s </title>\n</head>\n<body>\n<h2>Matlab Report Page : %s</h2>\n',FileName,FileName);
	
else
	fprintf(fPointer,'<html>\n<body>\n',FileName);
end
%now let's print figure
resol = ['-r' num2str(Resolution)];
print(FigHandle, '-djpeg70', [LongFigName '.jpg'],resol);
% and save figure
saveas(FigHandle, [LongFigName '.fig'], 'fig');

% now generate link to image 
fprintf(fPointer, '<img src="%s" alt="Figure %d">\n<br><br>\n',[FileName '/' FigName '.jpg'],FigIndex);

if isempty(Comment)
% obtain comment
prompt = ['Enter comments for figure' num2str(FigIndex)];

Text = inputdlg(prompt, 'Commnents dialog', 1, {['Figure ' num2str(FigIndex)]});
Comment = Text{1};
end
   
%puty comment in the html file
fprintf(fPointer,'<p>%s</p>\n<br>\n<hr>\n<br>\n',Comment);


fprintf(fPointer,'</body>\n</html>\n\n');
fclose(fPointer);

% fPointer=fopen([UserPath '/mrep/index.html'],'at+');
% fprintf(
%system('mrepindex');

if Preview
	web(LongFileName,'-browser');
end
