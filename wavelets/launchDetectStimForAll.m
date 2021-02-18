function launchDetectStimForAll(trucdeouf)

%trucdeouf contient : numéro de la souris, le premier et dernier recording à processer de cette souris,
%le numéro du stim channel pour ces recordings



for i=1:size(trucdeouf,1)
    foldername = num2str(trucdeouf(i,1));
    cd(foldername)
    
    for ii=trucdeouf(i,2):trucdeouf(i,3)
        
        if ii < 10
            FileBase = ['MVm' foldername '-0' num2str(ii)];
        else
            FileBase = ['MVm' foldername '-' num2str(ii)];
        end
             
        DetectStimOnMV(FileBase,trucdeouf(i,4),0.1,10,[],1);
        
    end
    cd ..
end
            