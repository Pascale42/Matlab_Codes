%function ColorRange = PlotCSD(y, t, CSDCh, step,  ColorRange, scale, SeeOutTrace, WhichCSDCh, csdtype, orient)
function ColorRange = PlotCSD(y,t, CSDCh, varargin)
if nargin<3
    CSDCh = [];
end
[step, ColorRange, scale, SeeOutTrace, WhichCSDCh,csdtype,orient] =DefaultArgs(varargin, {2,[], 1,2, [],'c',[]});

% 
% if iscell(CSDCh) & ~iscell(y)
%     ytmp = {};
%     ny = size(CSDCh,1);
%     nx = size(CSDCh,2);
%     for n1=1:ny
%         for n2 = 1:nx
%             ytmp{n1,n2} = y(:,CSDCh{n1,n2});
%         end
%     end
%     y = ytmp;
% end
% 
% 
% 
% if iscell(y) & length(y)==1
%     y = y{1};
% end
% 

if iscell(y)
    % recursively do each plot
    ny = size(y,1);
    nx = size(y,2);
    
    if isempty(ColorRange)
        ColorRange =-inf;
        for n1=1:ny
            for n2 = 1:nx
                chnum = size(y{n1,n2},2);
                ch=[step+1:chnum-step];
                %                y{n1,n2} = y{n1,n2} - repmat(mean(y{n1,n2}),size(y{n1,n2},1),1);
                
                csd = y{n1,n2}(:,ch+step) - 2*y{n1,n2}(:,ch) + y{n1,n2}(:,ch-step);
                csd = -csd/(step^2);
                ColorRange = max([abs(csd(:));ColorRange]);
            end;
        end
        ColorRange = [-ColorRange ColorRange];
    end
    
    for n1=1:ny
        for n2 = 1:nx
            x = y{n1,n2};
            if iscell(t)
                tcur = t{n1,n2};
            else
                tcur =t;
            end
            if ~isempty(CSDCh)
                csdch = CSDCh{n1,n2};
            else
                csdch =[];
            end
            
            if ~isempty(WhichCSDCh)
                if iscell(WhichCSDCh)
                    wcsdch = WhichCSDCh{n1,n2};
                else
                    wcsdch = WhichCSDCh;
                end
                
            else
                wcsdch =[];
            end
            
            
            %here just complicated selection of plot orientation from
            %orient or the data
            if isempty(orient) | (nx>1 & ny>1)
                pnx = nx; pny = ny; pn1 = n1; pn2 =n2;
            elseif strcmp(orient,'v')
                pnx = 1;
                if nx>ny
                    pny = nx; pn1=n1; pn2=n2;
                else
                    pny = ny; pn1=n2; pn2=n1;
                end
                
            elseif strcmp(orient,'h')
                pny = 1;
                if nx>ny
                    pnx = nx; pn1=n1; pn2=n2;
                else
                    pnx = ny; pn1=n2; pn2=n1;
                end
                
            end
            
            subplot(pny,pnx, pnx*(pn1-1)+pn2);
            
            PlotCSD(y{n1,n2}, tcur, csdch, step, ColorRange, scale, SeeOutTrace, wcsdch, csdtype);
        end
    end
    
else
    % single plot routine
    
    if ~isempty(CSDCh)
        MissingCh = setdiff([CSDCh(1):CSDCh(end)],CSDCh);
        y = FixEegChannels(y',MissingCh,'cubic');
        y = y';
    end
    if ~isempty(WhichCSDCh)
        y = y(:, WhichCSDCh);
    end
    if length(t)==2
        nt = size(y,1);
        win = 1000*(nt-1)/2/1250;
        tt = linspace(-win, win, nt);
        ti = find(tt>t(1) & tt<t(2));
        y = y(ti,:);
        t = tt(ti);
    end
    %    figure
    if isempty(ColorRange)
        chnum = size(y,2);
        ch=[step+1:chnum-step];
        csd = y(:,ch+step) - 2*y(:,ch) + y(:,ch-step);
        csd = -csd/(step^2);
        ColorRange = max(abs(csd(:)));
        ColorRange = [-ColorRange ColorRange];
    end
    
    CurSrcDns(y, t, csdtype, [], [], [], step, ColorRange);
    %axes
    if nargout>1
        csdout =  CurSrcDns(y, t, csdtype, [], [], [], step, ColorRange);
    end
    
    if SeeOutTrace==1
        curAxes = get(gca,'Position');
        curAxes = [curAxes(1) curAxes(2)+0.05*scale curAxes(3) curAxes(4)-0.1*scale];
        set(gca,'Position',curAxes);
    end
    
    
    
    
    colorbar
    if strcmp(csdtype,'c')
        hold on
        if SeeOutTrace
            nCh = size(y,2);
            AxesIncr = 2*step/(nCh-2*step);
            subplot 121
             CurSrcDns(y, t, csdtype, [], [], [], step, ColorRange);
            subplot 122; PlotManyCh(y,t, 1250, scale , 'k', 0, AxesIncr); % !!! PROBLEM HERE
        elseif SeeOutTrace==2
            PlotManyCh(y(:, step+1:end-step),t, 1250, scale , 'k', 0, 0);
        end
    end
end
