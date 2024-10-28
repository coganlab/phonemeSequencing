global TASK_DIR

global DUKEDIR
global BOX_DIR
%BOX_DIR='C:\Users\gcoga\Box';
%BOX_DIR='H:\Box Sync';
%RECONDIR='C:\Users\gcoga\Box\ECoG_Recon';


fDown = 100; %Downsampled Sampling Frequency
timeExtract = [-1 1.5];

Task=[];

Task.Name='Phoneme_Sequencing';
Subject = popTaskSubjectData(Task);

Task.Base.Name='Start';
Task.Base.Epoch='Start';
Task.Base.Time=[-0.5 0];


TASK_DIR=([BOX_DIR '/CoganLab/D_Data/' Task.Name]);
DUKEDIR=TASK_DIR;

Task.Conds(1).Name='All';
Task.Conds(1).Field(1).Name='AuditorywDelay';
Task.Conds(1).Field(1).Epoch='Auditory';
Task.Conds(1).Field(1).Time=[0.1 0.6];
Task.Conds(1).Field(2).Name='DelaywGo';
Task.Conds(1).Field(2).Epoch='Go';
Task.Conds(1).Field(2).Time=[-0.5 0];
Task.Conds(1).Field(3).Name='Response';
Task.Conds(1).Field(3).Epoch='ResponseStart';
Task.Conds(1).Field(3).Time=[-0.25 0.25];




%%
subjNames = extractfield(Subject,'Name');
% subjList = {'D76','D80','D81','D84','D85','D86','D88','D91','D92','D93','D94'...
%     'D95','D96','D97','D100','D102','D103'};
subjList = {'D102','D103'};
warning off;
for iSubject=1:length(subjList)
%    if exist([DUKEDIR '\Stats\timeGroup\' Subject(iSubject).Name '_' Task.Name '_' ...
%                  Task.Conds(1).Name '_' Task.Conds(1).Field(2).Name '_' Task.Base.Name '.mat'])
%        disp(['Skipping ' Subject(iSubject).Name])
%        continue
%    end
   subjId = find(ismember(subjNames,subjList{iSubject}));
   disp(['Processing ' Subject(subjId).Name])
    Trials=Subject(subjId).Trials;
    counterN=0;
    counterNR=0;
    noiseIdx=0;
    noResponseIdx=0;
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
   
    condIdx=ones(length(Subject(subjId).Trials),1);
    chanIdx=Subject(subjId).goodChannels;
    channelName = {Subject(subjId).ChannelInfo.Name};
    channelName = channelName(chanIdx);
    baseEpoch=Task.Base.Epoch;
    baseTimeRange=Task.Base.Time;    
    Trials=Subject(subjId).Trials(setdiff(1:length(Subject(subjId).Trials),noiseIdx));

    ieegBase=trialIEEG(Trials,chanIdx,baseEpoch,timeExtract.*1000);
    ieegBase = permute(ieegBase,[2,1,3]);
    fs = Subject(iSubject).Experiment.processing.ieeg.sample_rate;   
    ieegBaseStruct = ieegStructClass(ieegBase, fs, timeExtract, [1 fs/2], baseEpoch);

    clear ieegBase;
    ieegBaseCAR = extractCar(ieegBaseStruct);
    %ieegBaseCAR=extractCableCar(ieegBaseStruct,channelName);
    %ieegBaseCAR = ieegBaseStruct;
    clear ieegBaseStruct;
    [ieegBaseHG,ieegbasePower] = extractHiGamma(ieegBaseCAR,fDown,baseTimeRange);
    %normFactor = extractHGnormFactor(ieegBaseHG);
    if isfield(Trials,'ResponseStart')
        numField = 3;
    else
        numField = 2;
    end
    for iC=1:length(Task.Conds)
        for iF=1:numField
          %  if iC<=2
                Trials=Subject(subjId).Trials(setdiff(find(condIdx==iC),cat(2,noiseIdx,noResponseIdx)));
         %   else
         %       Trials=Subject(SN).Trials(setdiff(find(condIdx==iC),noiseIdx));
         %   end
            
            Epoch=Task.Conds(iC).Field(iF).Epoch;
            fieldTimeRange=Task.Conds(iC).Field(iF).Time;          
                        
            ieegField=trialIEEG(Trials,chanIdx,Epoch,timeExtract.*1000);
            ieegField = permute(ieegField,[2,1,3]);
            ieegFieldStruct = ieegStructClass(ieegField, fs, timeExtract, [1 fs/2], Epoch);
            
            clear ieegField
            % Common average referencing
            ieegFieldCAR = extractCar(ieegFieldStruct);
           % ieegFieldCAR = extractCableCar(ieegFieldStruct,channelName);
           %ieegFieldCAR = ieegFieldStruct;
            clear ieegFieldStruct;
            % High gamma extraction
            [ieegFieldHG,ieegFieldPower] = extractHiGamma(ieegFieldCAR,fDown,fieldTimeRange);    
            % Power comparison
            pChan = [];
            parfor iChan = 1:length(chanIdx)
                iChan
                pChan(iChan) = permtest(ieegFieldPower(iChan,:),ieegbasePower(iChan,:),10000);   
            end
            channelNames = {Subject(subjId).ChannelInfo(chanIdx).Name};
            
            if ~exist([DUKEDIR '\Stats\timeGroup\'])
                mkdir([DUKEDIR '\Stats\timeGroup\'])
            end
            save([DUKEDIR '\Stats\timeGroup\' Subject(subjId).Name '_' Task.Name '_' ...
                Task.Conds(iC).Name '_' Task.Conds(iC).Field(iF).Name '_' Task.Base.Name '.mat'],'pChan','channelNames','ieegFieldHG','ieegBaseHG');
         end
    end
end
