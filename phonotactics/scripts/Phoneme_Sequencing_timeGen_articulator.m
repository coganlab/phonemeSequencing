% global TASK_DIR
% 
% global DUKEDIR
 global BOX_DIR
%BOX_DIR='C:\Users\gcoga\Box';
% BOX_DIR='H:\Box Sync';
% RECONDIR='C:\Users\gcoga\Box\ECoG_Recon';
saveFolder = '\TempDecode\Articulator\sensorimotor_delay_sig_v2\';
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
Task.Conds(1).Name='SMC';
Task.Conds(1).Field(1).Name='AuditorywDelay';
Task.Conds(1).Field(1).Epoch='Auditory';
Task.Conds(1).Field(1).Time=[0 1.5];
Task.Conds(1).Field(2).Name='DelaywGo';
Task.Conds(1).Field(2).Epoch='Go';
Task.Conds(1).Field(2).Time=[-1 1.5];
Task.Conds(1).Field(3).Name='Response';
Task.Conds(1).Field(3).Epoch='ResponseStart';
Task.Conds(1).Field(3).Time=[-1 1];


% Decoder variable
numFold = 10;
varExplained = 80; 
% varExplained = [10:10:90]; Picking desired variance through
% cross-validation; Computationally expensive
numIter = 1; % Number of decoder iterations for stable results
timeWin = 0.2; % Window length for temporal generalization (s)
timeRes = 0.025; % Window hop for temporal generalization (s)

SNList=1:length(Subject);
%SNList=5;
 %% 
for iSubject=1:length(Subject)
  
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
    anatChan = (sensorymotorChan);
    disp(['Number of anatomical channels : ' num2str(sum(anatChan))]);
    
    % prodDelayElecs is loaded from the mat file in 'forKumar'
    % Come up with automated way to parse necessary channels
    channelName = {Subject(iSubject).ChannelInfo.Name};
    channelName(cellfun(@isempty,channelName)) = {'dummy'};
    channelName = channelName(chanIdx);
    [delayChan2select] = ismember(channelName,prodDelayElecs); 
    disp(['Number of delay activated channels : ' num2str(sum(delayChan2select))]);
    
    totalChan2Select = anatChan&delayChan2select;
    if(isempty(find(totalChan2Select)))
        disp('No significant anatomical channels with delay activation found; Iterating next subject');
        continue;
        % Forces the iteration for next subject;
    else
        disp(['Number of delay activated channels at selected anatomy: ' num2str(sum(totalChan2Select))]);
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
    ieegBaseCAR.data = ieegBaseCAR.data(totalChan2Select,:,:);
    ieegBaseHG = extractHiGamma(ieegBaseCAR,fDown,baseTimeRange);
    normFactorBase = extractHGnormFactor(ieegBaseHG);
    for iCond=1:length(Task.Conds)
        trials2Select = setdiff(find(condIdx==iCond),cat(2,noiseIdx,noResponseIdx));
        Trials=Subject(iSubject).Trials(trials2Select);
        trialInfo = Subject(iSubject).trialInfo(trials2Select);
        % Phoneme Sequence trial parsing
        phonemeTrial = phonemeSequenceTrialParser(trialInfo);
        
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
            % Normalized High gamma extraction
            ieegFieldCAR.data = ieegFieldCAR.data(totalChan2Select,:,:);
            ieegFieldHG = extractHiGamma(ieegFieldCAR,fDown,fieldTimeRange,normFactorBase,normType);
            
            
            % Initialize Decoder object
            phonDecode = decoderClass(numFold,varExplained,numIter);
            disp('1D temporal generalization');
            % 1D temporalGeneralization            
            decodeTimeStruct1D = phonDecode.tempGenClassify1D(ieegFieldHG,phonemeTrial.phonemeClass(:,1)',timeRes,timeWin,1:sum(totalChan2Select));
            % 2D temporalGeneralization            
            disp('2D temporal generalization');
            decodeTimeStruct2D = phonDecode.tempGenClassify2D(ieegFieldHG,phonemeTrial.phonemeClass(:,1)',timeRes,timeWin,1:sum(totalChan2Select),1:length(phonemeTrial.phonemeClass(:,1)'));
            
            channelsUsed = channelName(anatChan);
            % 
            disp('Saving..');
            if ~exist([DUKEDIR saveFolder]) %#ok<EXIST> 
                mkdir([DUKEDIR saveFolder])
            end
            save(fullfile([DUKEDIR saveFolder Subject(iSubject).Name '_' Task.Name '_' ...
                Task.Conds(iCond).Name '_' Task.Conds(iCond).Field(iField).Name '_' Task.Base.Name '_articulator_decoded.mat']),'decodeTimeStruct1D','decodeTimeStruct2D','channelsUsed');
        end
    end
end
