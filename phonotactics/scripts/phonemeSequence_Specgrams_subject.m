 global BOX_DIR

saveFolder = '\TempDecode\temporal_nosig\';
BOX_DIR = 'C:\Users\sd355\Box'
fDown = 200; %Downsampled Sampling Frequency
timeExtract = [-2.5 2.5];

Task=[];

Task.Name='Phoneme_Sequencing';
Subject = popTaskSubjectData(Task);

Task.Base.Name='Start';
Task.Base.Epoch='Start';
Task.Base.Time=[-0.5 0];


TASK_DIR=([BOX_DIR '/CoganLab/D_Data/' Task.Name]);
DUKEDIR=TASK_DIR;

% Task conditions
Task.Conds(1).Name='All';
Task.Conds(1).Field(1).Name='AuditorywDelay';
Task.Conds(1).Field(1).Epoch='Auditory';
Task.Conds(1).Field(1).Time=[-0.5 2];
% Task.Conds(1).Field(2).Name='DelaywGo';
% Task.Conds(1).Field(2).Epoch='Go';
% Task.Conds(1).Field(2).Time=[-0.5 1.5];
Task.Conds(1).Field(2).Name='Response';
Task.Conds(1).Field(2).Epoch='ResponseStart';
Task.Conds(1).Field(2).Time=[-1 1.5];
%%
for iSubject=31
     
    Subject(iSubject).Name
    Trials=Subject(iSubject).Trials;
    counterN=0;
    counterNR=0;
    noiseIdx=0;
    noResponseIdx=0;
    
    chanIdx=Subject(iSubject).goodChannels;

    anatName = {Subject(iSubject).ChannelInfo.Location};
    anatName(cellfun(@isempty,anatName)) = {'dummy'};
    anatName = anatName(chanIdx);
    sensorymotorChan = contains(anatName,'central');
    whiteChan = contains(anatName,'white');
    ifgChan = contains(anatName,'opercularis');
    frontalChan = contains(anatName,'front');
    temporalChan = contains(anatName,'temporal');
    anatChan = sensorymotorChan;
    disp(['Number of anatomical channels : ' num2str(sum(anatChan))]);
    
    % prodDelayElecs is loaded from the mat file in 'forKumar'
    % Come up with automated way to parse necessary channels
    channelName = {Subject(iSubject).ChannelInfo.Name};
    channelName(cellfun(@isempty,channelName)) = {'dummy'};
    channelName = channelName(chanIdx);
    [~,prodChan] = intersect(channelName,elecNameProductionClean); 
   
    chan2select = intersect(find(anatChan),prodChan);

     if(isempty(chan2select))
        disp('No channels found; Iterating next subject');
        continue;
        % Forces the iteration for next subject;
    end
    chanNameSubject = channelName(chan2select);
    
    for iTrials=1:length(Trials)
        if Trials(iTrials).Noisy==1
            noiseIdx(counterN+1)=iTrials;
            counterN=counterN+1;
        end
        if Trials(iTrials).NoResponse==1
            noResponseIdx(counterNR+1)=iTrials;
            counterNR=counterNR+1;
        end
    end
   
    condIdx=ones(length(Subject(iSubject).Trials),1);
    
    baseEpoch=Task.Base.Epoch;
    baseTimeRange=Task.Base.Time;    
    Trials=Subject(iSubject).Trials(setdiff(1:length(Subject(iSubject).Trials),noiseIdx));
    
    ieegBase=trialIEEG(Trials,chanIdx,baseEpoch,timeExtract.*1000);
    ieegBase = permute(ieegBase,[2,1,3]);
    fs = Subject(iSubject).Experiment.processing.ieeg.sample_rate;   
    ieegBaseStruct = ieegStructClass(ieegBase, fs, timeExtract, [1 fs/2], baseEpoch);
    clear ieegBase;
    ieegBaseCAR=extractCar(ieegBaseStruct);
    clear ieegBaseStruct;
    [specCarBase]= getSpectrograms(ieegBaseCAR.data(chan2select,:,:),[],ieegBaseCAR.tw,baseTimeRange,[1 500],[0 1],[0 1],[70 150],fs,0);
    meanFreqChanOut = extractSpecNorm(specCarBase,baseTimeRange,baseTimeRange);
    
    for iCond=1:length(Task.Conds)
        trials2Select = setdiff(find(condIdx==iCond),cat(2,noiseIdx,noResponseIdx));
        % Iterting through field: 'Auditory', 'Delay', 'Response'
        for iField=1:length(Task.Conds(iCond).Field)
         
            Epoch=Task.Conds(iCond).Field(iField).Epoch;
            fieldTimeRange=Task.Conds(iCond).Field(iField).Time;          
                        
            ieegField=trialIEEG(Trials,chanIdx,Epoch,timeExtract.*1000);
            ieegField = permute(ieegField,[2,1,3]);
            ieegFieldStruct = ieegStructClass(ieegField, fs, timeExtract, [1 fs/2], Epoch);
            clear ieegField
            % Common average referencing
            ieegFieldCAR = extractCar(ieegFieldStruct);
            clear ieegFieldStruct;
            [specCarField]= getSpectrograms(ieegFieldCAR.data(chan2select,:,:),[],ieegFieldCAR.tw,fieldTimeRange,[1 500],[0 1],[0 1],[70 150],fs,0);
            for iChan = 1:length(chan2select)
                figure;
                specChanMap(specCarField,[],1:length(chan2select),[],fieldTimeRange,[-1.5 -1],[-0.25 0.25],[1 500],[70 150],[-2 2],iChan,meanFreqChanOut);
                xlabel(['Time from ' Epoch ' onset (s)'])
                ylabel('Frequency (Hz)')
                title(chanNameSubject(iChan))
            end
        
        end
        
    end
end