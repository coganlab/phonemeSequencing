 global BOX_DIR
%BOX_DIR='C:\Users\gcoga\Box';
% BOX_DIR='H:\Box Sync';
% RECONDIR='C:\Users\gcoga\Box\ECoG_Recon';
saveFolder = '\TempDecode\temporal_nosig\';
BOX_DIR = 'C:\Users\sd355\Box'
fDown = 200; %Downsampled Sampling Frequency
timeExtract = [-2 2];

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
for iSubject=6
     
    Subject(iSubject).Name
    Trials=Subject(iSubject).Trials;
    counterN=0;
    counterNR=0;
    noiseIdx=0;
    noResponseIdx=0;
    
    %chanIdx=Subject(iSubject).goodChannels;
    chanIdx = 1:length(Subject(iSubject).ChannelInfo);
    anatName = {Subject(iSubject).ChannelInfo.Location};
    anatName(cellfun(@isempty,anatName)) = {'dummy'};
    anatName = anatName(chanIdx);
    sensorymotorChan = contains(anatName,'central');
    whiteChan = contains(anatName,'white');
    ifgChan = contains(anatName,'opercularis');
    frontalChan = contains(anatName,'front');
    temporalChan = contains(anatName,'temporal');
    anatChan = ones(1,length(anatName));
    disp(['Number of anatomical channels : ' num2str(length(anatChan))]);
    
    % prodDelayElecs is loaded from the mat file in 'forKumar'
    % Come up with automated way to parse necessary channels
    
    channelName = {Subject(iSubject).ChannelInfo.Name};
    channelName(cellfun(@isempty,channelName)) = {'dummy'};
    channelName = channelName(chanIdx);
    %muscleChannelsName = channelName(muscleChannels);
    [~,delayChan] = intersect(channelName,elecNameProductionClean); 
   
    chan2select = intersect(find(anatChan),delayChan);
    %chan2select = muscleChannels;
     if(isempty(chan2select))
        disp('No channels found; Iterating next subject');
        continue;
        % Forces the iteration for next subject;
    end
    chanNameSubject = channelName(chan2select);

%%
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
    [NumTrials,goodtrials] = remove_bad_trials(ieegBaseCAR.data,10);
    goodTrialsCommon = extractCommonTrials(goodtrials);
    waveSpecBase = getWaveletScalogram(ieegBaseCAR.data(chan2select,goodTrialsCommon,:),ieegBaseCAR.fs);
    meanFreqChanOut = extractSpecNorm(waveSpecBase.spec,timeExtract,baseTimeRange);
    
    for iCond=length(Task.Conds)
        trials2Select = setdiff(find(condIdx==iCond),cat(2,noiseIdx,noResponseIdx));
        % Iterting through field: 'Auditory', 'Delay', 'Response'
        for iField=1:length(Task.Conds(iCond).Field)
            Trials = Subject(iSubject).Trials(trials2Select);
            Epoch=Task.Conds(iCond).Field(iField).Epoch;
            fieldTimeRange=Task.Conds(iCond).Field(iField).Time;          
                        
            ieegField=trialIEEG(Trials,chanIdx,Epoch,timeExtract.*1000);
            ieegField = permute(ieegField,[2,1,3]);
            ieegFieldStruct = ieegStructClass(ieegField, fs, timeExtract, [1 fs/2], Epoch);
            clear ieegField
            % Common average referencing
            ieegFieldCAR = extractCar(ieegFieldStruct);
            clear ieegFieldStruct;
            [NumTrials,goodtrials] = remove_bad_trials(ieegFieldCAR.data,10);
            goodTrialsCommon = extractCommonTrials(goodtrials);
            waveSpecField = getWaveletScalogram(ieegFieldCAR.data(chan2select,goodTrialsCommon,:),ieegFieldCAR.fs);
            for iChan = 1:length(chan2select)
                figure;

                subplot(2,1,1);
                timeVect = linspace(ieegFieldCAR.tw(1),ieegFieldCAR.tw(2),size(ieegFieldCAR.data,3));
                
                plot(timeVect,squeeze(ieegFieldCAR.data(chan2select(iChan),goodTrialsCommon,:)),'color',[0 0 0] +0.75);
                hold on;
                plot(timeVect,mean(squeeze(ieegFieldCAR.data(chan2select(iChan),goodTrialsCommon,:)),1),'color',[0 0 0]);
                xlim(fieldTimeRange);
               

                subplot(2,1,2);
                specChanMap(waveSpecField.spec,[],1:length(chan2select),[],timeExtract,[-1.5 -1],[-0.25 0.25],[1 500],[70 150],[-2 2],iChan,meanFreqChanOut);
                set(gca,'YTick',1:4:length(waveSpecField.fscale));
                set(gca,'YTickLabels',waveSpecField.fscale(1:4:end));
                xlim(fieldTimeRange);
                xlabel(['Time from ' Epoch ' onset (s)'])
                ylabel('Frequency (Hz)')
                title(chanNameSubject(iChan))
            end
        
        end
        
    end
end
