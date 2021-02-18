%
% Pyramid-Interneuron units / Field Multitaper Coherence Density
%
% function OutArgs = UCTFCoherenceDensity(FileBase, State, fMode, Overwrite, SaveFig);
%
% !! COMPUTE UFCoherenceDensity BEFORE !!
% !! COMPUTE NeuroClass AND myCellType BEFORE !!
%
% State: THE SWS GAMTHE SWSTHE
% SaveFig: 1 to save report
% fMode: 'compute' and 'display' modes


function OutArgs = UCTFCoherenceDensity(FileBase, State, fMode, Overwrite, SaveFig);

switch fMode

    case 'compute'

        if exist([FileBase '.' mfilename State '.mat'])>0 & ~Overwrite
            sprintf('Already computed!');return
        end

        load([FileBase '.UFCoherenceDensity' State '.mat']);
        load([FileBase '.myCellType.mat']);

        % Creates a new MAP to sort Pyr form IN from nd;

        n=length(OutArgs.Map); % !! columns are indx-clu-shk !!
        m=length(MapCT); % !! columns are indx-shk-clu !!
        MAP=[];
        map=[];
        for i=1:n
            for j=1:m
                if OutArgs.Map(i,3) == MapCT(j,2) & OutArgs.Map(i,2) == MapCT(j,3);
                    map=[OutArgs.Map(i,:) MapCT(j,4)];
                end
            end
            MAP=[MAP;map];
        end
        MAP=[MAP(:,1) MAP(:,3) MAP(:,2) MAP(:,4)]; % put back order indx-shk-clu-CT

        OutArgsCT= OutArgs;
        OutArgsCT.MAP = MAP;
        OutArgsCT = rmfield(OutArgsCT, 'Map');
        OutArgsCT = rmfield(OutArgsCT, 'ufchd');

        nCh = length(OutArgs.Channels);

        % Now Reassign ufchd fields for Pyr IN and Nd;
        % y
        OutArgsCT.ufchd.Pyr.y = OutArgs.ufchd.y(:,[1:nCh],nCh+(find(MAP(:,4) == 1)) );
        OutArgsCT.ufchd.IN.y = OutArgs.ufchd.y(:,[1:nCh],nCh+(find(MAP(:,4) == -1)) );
        OutArgsCT.ufchd.Nd.y = OutArgs.ufchd.y(:,[1:nCh],nCh+(find(MAP(:,4) == 0)) );

        % f
        OutArgsCT.ufchd.f = OutArgs.ufchd.f;

        % phi
        OutArgsCT.ufchd.Pyr.phi = OutArgs.ufchd.phi(:,[1:nCh],nCh+(find(MAP(:,4) == 1)) );
        OutArgsCT.ufchd.IN.phi = OutArgs.ufchd.phi(:,[1:nCh],nCh+(find(MAP(:,4) == -1)) );
        OutArgsCT.ufchd.Nd.phi = OutArgs.ufchd.phi(:,[1:nCh],nCh+(find(MAP(:,4) == 0)) );

        % yerr
        OutArgsCT.ufchd.Pyr.yerr = OutArgs.ufchd.yerr(:,[1:nCh],nCh+(find(MAP(:,4) == 1)),: );
        OutArgsCT.ufchd.IN.yerr = OutArgs.ufchd.yerr(:,[1:nCh],nCh+(find(MAP(:,4) == -1)),: );
        OutArgsCT.ufchd.Nd.yerr = OutArgs.ufchd.yerr(:,[1:nCh],nCh+(find(MAP(:,4) == 0)),: );

        % phierr
        OutArgsCT.ufchd.Pyr.phierr = OutArgs.ufchd.phierr(:,[1:nCh],nCh+(find(MAP(:,4) == 1)) );
        OutArgsCT.ufchd.IN.phierr = OutArgs.ufchd.phierr(:,[1:nCh],nCh+(find(MAP(:,4) == -1)) );
        OutArgsCT.ufchd.Nd.phierr = OutArgs.ufchd.phierr(:,[1:nCh],nCh+(find(MAP(:,4) == 0)) );

        % phloc
        OutArgsCT.ufchd.Pyr.phloc = OutArgs.ufchd.phloc(:,[1:nCh],nCh+(find(MAP(:,4) == 1)) );
        OutArgsCT.ufchd.IN.phloc = OutArgs.ufchd.phloc(:,[1:nCh],nCh+(find(MAP(:,4) == -1)) );
        OutArgsCT.ufchd.Nd.phloc = OutArgs.ufchd.phloc(:,[1:nCh],nCh+(find(MAP(:,4) == 0)) );

        % pow
        OutArgsCT.ufchd.Pyr.phloc = OutArgs.ufchd.phloc(:,[1:nCh],nCh+(find(MAP(:,4) == 1)) );
        OutArgsCT.ufchd.IN.phloc = OutArgs.ufchd.phloc(:,[1:nCh],nCh+(find(MAP(:,4) == -1)) );
        OutArgsCT.ufchd.Nd.phloc = OutArgs.ufchd.phloc(:,[1:nCh],nCh+(find(MAP(:,4) == 10)) );

        save([FileBase '.' mfilename State '.mat'], 'OutArgsCT');

    case 'display'

        load([FileBase '.' mfilename State '.mat']);

        nChannels = length(OutArgsCT.Channels);

        if isempty(OutArgsCT.ufchd.Pyr.y) == 0
            figure(100)
            for i=1:nChannels
                subplot(2,2,i);
                imagesc(OutArgsCT.ufchd.f,[1:length(OutArgsCT.MAP(find(OutArgsCT.MAP(:,4) == 1)))],sq(OutArgsCT.ufchd.Pyr.y(:,i,:))');
                set(gca,'YTick',[1:length(OutArgsCT.MAP(find(OutArgsCT.MAP(:,4) == 1)))]);
                set(gca,'YTickLabel',num2str(OutArgsCT.MAP((find(OutArgsCT.MAP(:,4) == 1)),2:3)));
                xlabel('Hz'); title(['Units Pyr - Field (ch ' num2str(OutArgsCT.Channels(i)) ') coherence density for ' num2str(State)]);
            end
        end
        if SaveFig > 0;
            reportfig(gcf,['UCTFCoherenceDensity_' State], 0, ['File ' OutArgsCT.Par.FileName ', State ' State ', ' num2str(OutArgsCT.FreqRange) ' Hz'],150);
        end

        if isempty(OutArgsCT.ufchd.IN.y) == 0
            figure(101)
            for i=1:nChannels
                subplot(2,2,i);
                imagesc(OutArgsCT.ufchd.f,[1:length(OutArgsCT.MAP(find(OutArgsCT.MAP(:,4) == -1)))],sq(OutArgsCT.ufchd.IN.y(:,i,:))');
                set(gca,'YTick',[1:length(OutArgsCT.MAP(find(OutArgsCT.MAP(:,4) == -1)))]);
                set(gca,'YTickLabel',num2str(OutArgsCT.MAP((find(OutArgsCT.MAP(:,4) == -1)),2:3)));
                xlabel('Hz'); title(['Units IN - Field (ch ' num2str(OutArgsCT.Channels(i)) ') coherence density for ' num2str(State)]);
            end
        end
        if SaveFig > 0;
            reportfig(gcf,['UCTFCoherenceDensity_' State], 0, ['File ' OutArgsCT.Par.FileName ', State ' State ', ' num2str(OutArgsCT.FreqRange) ' Hz'],150);
        end
        if isempty(OutArgsCT.ufchd.Nd.y) == 0
            figure(102)
            for i=1:nChannels
                subplot(2,2,i);
                imagesc(OutArgsCT.ufchd.f,[1:length(OutArgsCT.MAP(find(OutArgsCT.MAP(:,4) == 0)))],sq(OutArgsCT.ufchd.Nd.y(:,i,:))');
                set(gca,'YTick',[1:length(OutArgsCT.MAP(find(OutArgsCT.MAP(:,4) == 0)))]);
                set(gca,'YTickLabel',num2str(OutArgsCT.MAP((find(OutArgsCT.MAP(:,4) == 0)),2:3)));
                xlabel('Hz'); title(['Units Nd - Field (ch ' num2str(OutArgsCT.Channels(i)) ') coherence density for ' num2str(State)]);
            end
        end

        if SaveFig > 0;
            reportfig(gcf,['UCTFCoherenceDensity_' State], 0, ['File ' OutArgsCT.Par.FileName ', State ' State ', ' num2str(OutArgsCT.FreqRange) ' Hz'],150);
        end
end