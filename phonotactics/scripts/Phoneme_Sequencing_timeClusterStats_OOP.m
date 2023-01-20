global TASK_DIR

global DUKEDIR
global BOX_DIR
%BOX_DIR='C:\Users\gcoga\Box';
%BOX_DIR='H:\Box Sync';
%RECONDIR='C:\Users\gcoga\Box\ECoG_Recon';


fDown = 100; %Downsampled Sampling Frequency
timeExtract = [-5 5];

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
Task.Conds(1).Field(1).Time=[-0.5 1.5];
Task.Conds(1).Field(2).Name='DelaywGo';
Task.Conds(1).Field(2).Epoch='Go';
Task.Conds(1).Field(2).Time=[-1 1.5];
Task.Conds(1).Field(3).Name='Response';
Task.Conds(1).Field(3).Epoch='ResponseStart';
Task.Conds(1).Field(3).Time=[-1 1];
Task.Conds(1).Field(4).Name='StartCue';
Task.Conds(1).Field(4).Epoch='Start';
Task.Conds(1).Field(4).Time=[-0.5 1];



SNList=1:length(Subject);
%SNList=7;
for iSN=1:length(SNList)
    SN=SNList(iSN);
    Trials=Subject(SN).Trials;
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
   
    condIdx=ones(length(Subject(SN).Trials),1);
    chanIdx=Subject(SN).goodChannels;
    baseEpoch=Task.Base.Epoch;
    baseTimeRange=Task.Base.Time;    
    Trials=Subject(SN).Trials(setdiff(1:length(Subject(SN).Trials),noiseIdx));

    ieegBase=trialIEEG(Trials,chanIdx,baseEpoch,timeExtract.*1000);
    ieegBase = permute(ieegBase,[2,1,3]);
    fs = Subject(SN).Experiment.processing.ieeg.sample_rate;   
    ieegBaseStruct = ieegStructClass(ieegBase, fs, timeExtract, [1 fs/2], baseEpoch);

    clear ieegBase;
    ieegBaseCAR=extractCar(ieegBaseStruct);
    clear ieegBaseStruct;
    ieegBaseHG = extractHiGamma(ieegBaseCAR,fDown,baseTimeRange);
   
    for iC=1:length(Task.Conds)
        for iF=1:length(Task.Conds(iC).Field)
          %  if iC<=2
                Trials=Subject(SN).Trials(setdiff(find(condIdx==iC),cat(2,noiseIdx,noResponseIdx)));
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
            clear ieegFieldStruct;
            % High gamma extraction
            ieegFieldHG = extractHiGamma(ieegFieldCAR,fDown,fieldTimeRange);    
            % Time Series Permutation cluster
            chanSig = extractTimePermCluster(ieegFieldHG,ieegBaseHG);         
            channelNames = {Subject(SN).ChannelInfo(chanIdx).Name};
            
%             if ~exist([DUKEDIR '\Stats\timePerm\'])
%                 mkdir([DUKEDIR '\Stats\timePerm\'])
%             end
%             save([DUKEDIR '\Stats\timePerm\' Subject(SN).Name '_' Task.Name '_' ...
%                 Task.Conds(iC).Name '_' Task.Conds(iC).Field(iF).Name '_' Task.Base.Name '.mat'],'chanSig','channelNames','ieegFieldHG','ieegBaseHG');
        end
    end
end
