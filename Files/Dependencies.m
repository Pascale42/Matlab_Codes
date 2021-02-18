function [] = Dependencies(fullPath)
  %% DEPENDENCIES ON CURRENT FOLDER
  %  Show dependencies between each files of the current folder
  %  fullPath = 'full' to have the full path of the dependence
  %  fullPath = '' to have only the dependence's name
  %  Change MatlabToolBoxDir to avoid thoses dependencies
  %  Mathieu.GOLOS@gmail.com
  %  23/01/2013


  %% Parameters
  MatlabToolBoxDir = '/Applications/MATLAB_R2012a/';
  listing = dir;
  expr = '[^\'']\''([^\''\n\r]+(?:\w\.(?:mat|txt)){1})\''[^\'']';


  %% Loop on directory's files
  for i=1:size(listing,1)
    if ~ isdir(listing(i).name)
      % Dependencies with *.m
      fcnName = listing(i).name;
      fcnList = depfun(fcnName, '-toponly', '-quiet');
      listInd = strmatch(MatlabToolBoxDir, fcnList);
      strList = char(fcnList(setdiff(1:numel(fcnList),listInd)));

      % Dependencies with *.mat and *.txt
      fid = fopen(fcnName);
      fcnText = fscanf(fid, '%c');
      fclose(fid);
      dataFiles = regexp(fcnText,expr,'tokens');
      dataFiles = unique([dataFiles{:}]).';

      % Display if there is at least one dependence
      if (size(strList, 1) > 1) || (size(dataFiles,1) > 0)
        disp(listing(i).name)
        % Dependencies with *.m
        for j=1 : size(strList, 1)
          idx = strfind(strList(j, :), '\');
          fname = strtrim(strList(j, idx(end) +1 : end));
          if ~ strcmp(strtrim(listing(i).name), fname)
            if fullPath
              disp(['    ', strList(j,:)])
            else
              disp(['    ', fname])
            end
          end
        end
        % Dependencies with *.mat and *.txt
        for j=1 : size(dataFiles, 1)
          disp(['    ', char(dataFiles(j))])
        end
      end
    end
  end
end

