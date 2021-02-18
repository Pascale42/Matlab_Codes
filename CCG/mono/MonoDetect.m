function mono =  MonoDetect(FileBase,varargin)
%function MonoDetect(FileBase,Overwrite,ElLoc,Control,Report,Periods)
mono=struct([]);
par = LoadPar([FileBase '.par']);
[Overwrite,ElLoc,Control,Report,Periods] = DefaultArgs(varargin,{1,'c',1,0, {{'REM','SWS','RUN'}}});

if isstr(ElLoc(1))
    El = find(strcmp(par.ElecLoc,ElLoc));
else
    El = ElLoc;
end

if isempty(El)
    fprintf('No cells on electrodes from %s\n',ElLoc);
    return;
end
Fs = 1e6/par.SampleTime;
s2sFn = [FileBase '.s2s'];
if FileExists(s2sFn)
    load(s2sFn,'-MAT');
else
    s2s = Spikes2Spikes(FileBase,1,Periods,El);
end

SaveFn = [FileBase '.mono'];

if ~FileExists(SaveFn) | Overwrite

    myClu = find(ismember(s2s.ElClu(:,1),El)); %cells from your location
    nClu = length(myClu);
    %for where = 1:length(Periods) %loop through states - later
    % stored in s2s file
    %             s2s(where).cID = cID(myClu);
    %                 s2s(where).ElClu = ClustByEl(myClu,:);
    %                 s2s(where).ccg = ccg;
    %                 s2s(where).tbin =t;
    %                 s2s(where).State = Ext{where};

    %now go through cells of your location
    cnt=1;
    plotscnt=0;
    for i=1:nClu
        c1 = myClu(i);
        for j=i+1:nClu

            c2 = myClu(j);
       
            xc = sq(s2s.ccg(:,c1,c2,4)); % look at total counts for now
            middle = (length(xc)-1)/2+1;
            if s2s.ElClu(myClu(i),1)==s2s.ElClu(c2,1)
                sameel =1;
            else
                sameel =0;
            end

            if isempty(xc)
                continue;
            end
            [l,r, thr] = TestMono(xc,sameel,Control);
            if (l==2 && r==1) || (l==1 && r==2)
                mutual =1;
                
            else
                mutual =0;
            end
            if l>0

                mono(cnt).from = [s2s.cID(c2) s2s.ElClu(c2,:)];
                mono(cnt).to = [s2s.cID(c1) s2s.ElClu(c1,:)];
                mono(cnt).type = l;
                mono(cnt).ccg = xc(middle+[-50:50])';
                mono(cnt).thr = thr;
                mono(cnt).mutual = mutual;
                cnt=cnt+1;
                
            end
            if r>0
                mono(cnt).from = [s2s.cID(c1) s2s.ElClu(c1,:)];
                mono(cnt).to = [s2s.cID(c2) s2s.ElClu(c2,:) ];
                mono(cnt).type = r;
                mono(cnt).ccg = xc(middle+[-50:50])';
                mono(cnt).thr = thr;
                mono(cnt).mutual = mutual;
                cnt=cnt+1;
            end
            if l>0 || r>0
                plotscnt=plotscnt+1;
                if Control
					
	        	    pairstr = [FileBase '   :   ' num2str(s2s.ElClu(c1,:)) ' -> ' num2str(s2s.ElClu(c1,:))];
    		        title([pairstr ':  ' num2str([l r])]);

                    waitforbuttonpress

                    sel = get(gcf,'SelectionType');
                    fprintf('Left mouse click - go further, right click - keyboard command\n');

                    switch sel
                        case 'normal'
                            continue
                        case 'alt'
                            keyboard
                    end
                end


            end

        end
    end
    nPairs = cnt-1;
    if nPairs==0
            fprintf('no monopairs in this file\n');
            return;
    end
    
    
    if Report
        
        figure(123)
        clf
        mutFlag =0;
        cnt=1;
        for i=1:nPairs
            if ~mutFlag
                %first time in the mutual pair
                subplotfit(cnt,plotscnt);
                bar([-50:50],mono(i).ccg);axis tight
                Lines([],mono(i).thr,'r');
                from=num2str(mono(i).from(2:end));
                 to = num2str(mono(i).to(2:end));
                 if mono(i).mutual==1
                     dirstr = 'bidirect';
                 else
                     dirstr = 'unidirect';
                 end
                titstr =sprintf( '%s -> %s, %s',from, to, dirstr);
                title(titstr);
                cnt=cnt+1;
                if mono(1).mutual==1
                    mutFlag=1;
                end
            else
                mutFlag=0;
            end
        end
        
        reportfig(gcf,['mono_' ElLoc],0,FileBase,100);
    end

end




