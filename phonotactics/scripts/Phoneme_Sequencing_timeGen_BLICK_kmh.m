% global TASK_DIR
% 
% global DUKEDIR
% global BOX_DIR
%BOX_DIR='C:\Users\gcoga\Box';
% BOX_DIR='H:\Box Sync';
% RECONDIR='C:\Users\gcoga\Box\ECoG_Recon';


fDown = 200; %Downsampled Sampling Frequency
timeExtract = [-2 2];

Task=[];

Task.Name='Phoneme_Sequencing';
Subject = popTaskSubjectData(Task);

Task.Base.Name='Start';
Task.Base.Epoch='Start';
Task.Base.Time=[-0.5 0];

BOX_DIR = 'C:\Users\kmh156\Box';
TASK_DIR=([BOX_DIR '/CoganLab/D_Data/' Task.Name]);
DUKEDIR=TASK_DIR;

% Task conditions
Task.Conds(1).Name='All';
Task.Conds(1).Field(1).Name='AuditorywDelay';
Task.Conds(1).Field(1).Epoch='Auditory';
Task.Conds(1).Field(1).Time=[0 1.5];
Task.Conds(1).Field(2).Name='DelaywGo';
Task.Conds(1).Field(2).Epoch='Go';
Task.Conds(1).Field(2).Time=[-1 1.5];
Task.Conds(1).Field(3).Name='Response';
Task.Conds(1).Field(3).Epoch='ResponseStart';
Task.Conds(1).Field(3).Time=[-1 1];

% Regression variable
% Use one of the following: 'POS1', 'POS1', etc.
phonRegress = 'BLICK';
loop_count = 0;  % Kassie generated variable, =38 when error occurs
% Decoder variable
numFold = 10;
varExplained = 80; 
% varExplained = [10:10:90]; Picking desired variance through
% cross-validation; Computationally expensive
timeWin = 0.25; % Window length for temporal generalization (s)
timeRes = 0.05; % Window hop for temporal generalization (s)

SNList=1:length(Subject);
%SNList=5;
for iSN=1:length(SNList)
    SN=SNList(iSN);
    Subject(SN).Name
    Trials=Subject(SN).Trials;
    counterN=0;
    counterNR=0;
    noiseIdx=0;
    noResponseIdx=0;
    
    chanIdx=Subject(SN).goodChannels;
    
    % prodDelayElecs is loaded from the mat file in 'forKumar'
    % Come up with automated way to parse necessary channels
    loop_count = loop_count + 1;
    channelNames = {Subject(SN).ChannelInfo.Name};
    channelNames(cellfun(@isempty,channelNames)) = {'dummy'};
    channelNames = channelNames(chanIdx);
    [~,chan2select] = intersect(channelNames,prodDelayElecs); 
    if(isempty(chan2select))
        disp('No channels found; Iterating next subject');
        continue;
        % Forces the iteration for next subject;
    end
    
    
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
   
    condIdx=ones(length(Subject(SN).Trials),1);
    
    baseEpoch=Task.Base.Epoch;
    baseTimeRange=Task.Base.Time;    
    Trials=Subject(SN).Trials(setdiff(1:length(Subject(SN).Trials),noiseIdx));
    
    
    
    
    
    ieegBase=trialIEEG(Trials,chanIdx,baseEpoch,timeExtract.*1000);
    ieegBase = permute(ieegBase,[2,1,3]);
    fs = Subject(SN).Experiment.recording.sample_rate;   
    ieegBaseStruct = ieegStructClass(ieegBase, fs, timeExtract, [1 fs/2], baseEpoch);
    clear ieegBase;
    ieegBaseCAR=extractCar(ieegBaseStruct);
    clear ieegBaseStruct;
    ieegBaseHG = extractHiGamma(ieegBaseCAR,fDown,baseTimeRange);
    normFactorBase = extractHGnormFactor(ieegBaseHG);
    for iC=1:length(Task.Conds)
        trials2Select = setdiff(find(condIdx==iC),cat(2,noiseIdx,noResponseIdx));
        Trials=Subject(SN).Trials(trials2Select);
        trailInfo = Subject(SN).trialInfo(trials2Select);
        % Phoneme Sequence trial parsing
        phonemeTrial = phonemeSequenceTrialParser(trailInfo);
        
        % Iterting through field: 'Auditory', 'Delay', 'Production'
        for iF=1:length(Task.Conds(iC).Field)
         
            Epoch=Task.Conds(iC).Field(iF).Epoch;
            fieldTimeRange=Task.Conds(iC).Field(iF).Time;          
                        
            ieegField=trialIEEG(Trials,chanIdx,Epoch,timeExtract.*1000);
            ieegField = permute(ieegField,[2,1,3]);
            ieegFieldStruct = ieegStructClass(ieegField, fs, timeExtract, [1 fs/2], Epoch);
            clear ieegField
            % Common average referencing
            ieegFieldCAR = extractCar(ieegFieldStruct);
            clear ieegFieldStruct;
            % Normalized High gamma extraction
            ieegFieldHG = extractHiGamma(ieegFieldCAR,fDown,fieldTimeRange,normFactorBase);
            
            
            % Initialize Decoder object
            phonDecode = phonemeDecoderClass(numFold,varExplained);
            disp('1D temporal generalization');
            % 1D temporalGeneralization            
            decodeTimeStruct1D = tempGenRegress1D(phonDecode,ieegFieldHG,phonemeTrial,phonRegress,timeRes,timeWin,find(chan2select));
            % 2D temporalGeneralization            
            disp('2D temporal generalization');
            decodeTimeStruct2D = tempGenRegress2D(phonDecode,ieegFieldHG,phonemeTrial,phonRegress,timeRes,timeWin,find(chan2select));
            
            channelsUsed = channelNames(chan2select);
            % 
            disp('Saving..');
            if ~exist([DUKEDIR '\TempDecode\tempGen\'])
                mkdir([DUKEDIR '\TempDecode\tempGen\'])
            end
            save([DUKEDIR '\TempDecode\tempGen\' Subject(SN).Name '_' Task.Name '_' ...
                Task.Conds(iC).Name '_' Task.Conds(iC).Field(iF).Name '_' Task.Base.Name '_' phonRegress 'decoded.mat'],'decodeTimeStruct1D','decodeTimeStruct2D','channelsUsed');
        end
    end
end
