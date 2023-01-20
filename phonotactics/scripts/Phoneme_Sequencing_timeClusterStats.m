global TASK_DIR
global experiment
global DUKEDIR
global BOX_DIR
global RECONDIR
%BOX_DIR='C:\Users\gcoga\Box';
BOX_DIR='H:\Box Sync';
%RECONDIR='C:\Users\gcoga\Box\ECoG_Recon';
RECONDIR=[BOX_DIR '\ECoG_Recon'];

Task=[];

Task.Name='Phoneme_Sequencing';
Task.Base.Name='Start';
Task.Base.Epoch='Start';
Task.Base.Time=[-500 0];


TASK_DIR=([BOX_DIR '/CoganLab/D_Data/' Task.Name]);
DUKEDIR=TASK_DIR;

Task.Conds(1).Name='All';
Task.Conds(1).Field(1).Name='AuditorywDelay';
Task.Conds(1).Field(1).Epoch='Auditory';
Task.Conds(1).Field(1).Time=[-500 1500];
Task.Conds(1).Field(2).Name='DelaywGo';
Task.Conds(1).Field(2).Epoch='Go';
Task.Conds(1).Field(2).Time=[-1000 1500];
Task.Conds(1).Field(3).Name='Response';
Task.Conds(1).Field(3).Epoch='ResponseStart';
Task.Conds(1).Field(3).Time=[-1000 1000];
Task.Conds(1).Field(4).Name='StartCue';
Task.Conds(1).Field(4).Epoch='Start';
Task.Conds(1).Field(4).Time=[-500 1000];

Subject = popTaskSubjectData(Task);

SNList=1:length(Subject);
SNList=31:36;
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
    %  [condIdx noiseIdx noResponseIdx]=SentenceRepConds(Subject(SN).Trials);
    %chanIdx=setdiff(Subject(SN).goodChannels,Subject(SN).WM);
    condIdx=ones(length(Subject(SN).Trials),1);
    chanIdx=Subject(SN).goodChannels;
    baseEpoch=Task.Base.Epoch;
    baseTimeRange=Task.Base.Time;
    baseTimeRange(1)=baseTimeRange(1)-500;
    baseTimeRange(2)=baseTimeRange(2)+500;
    Trials=Subject(SN).Trials(setdiff(1:length(Subject(SN).Trials),noiseIdx));
    ieegBase=trialIEEG(Trials,chanIdx,baseEpoch,baseTimeRange);
    ieegBaseCAR=ieegBase-mean(ieegBase,2);
    TimeLengthBase=baseTimeRange(2)-baseTimeRange(1);
    for iC=1:length(Task.Conds)
        for iF=1:length(Task.Conds(iC).Field)
          %  if iC<=2
                Trials=Subject(SN).Trials(setdiff(find(condIdx==iC),cat(2,noiseIdx,noResponseIdx)));
         %   else
         %       Trials=Subject(SN).Trials(setdiff(find(condIdx==iC),noiseIdx));
         %   end
            
            Epoch=Task.Conds(iC).Field(iF).Epoch;
            TimeRange=Task.Conds(iC).Field(iF).Time;
            TimeRange(1)=TimeRange(1)-500;
            TimeRange(2)=TimeRange(2)+500;
            TimeLength=TimeRange(2)-TimeRange(1);
            
            ieeg=trialIEEG(Trials,chanIdx,Epoch,TimeRange);
            ieegCAR=ieeg-mean(ieeg,2);
            
            ieegCARHGZ=zeros(size(ieegCAR,2),size(ieegCAR,1),TimeLength./10);
            ieegBaseCARHGZ=zeros(size(ieegBaseCAR,2),size(ieegBaseCAR,1),TimeLengthBase./10);
            
            ieegCARHG=zeros(size(ieegCAR,2),size(ieegCAR,1),TimeLength./10);
            ieegBaseCARHG=zeros(size(ieegBaseCAR,2),size(ieegBaseCAR,1),TimeLengthBase./10);
            
            for iChan=1:size(ieegCAR,2)
                %   [~,dh2] = EcogExtractHighGammaTrial(sq(ieegCAR(:,iChan,:)), experiment.recording.sample_rate, 100, [70 120],[-1 2],[-1 2],[]);
                [~,hgsig] = EcogExtractHighGammaTrial(sq(ieegCAR(:,iChan,:)), Subject(SN).Experiment.recording.sample_rate, 100, [70 150],TimeRange,TimeRange,[]);
                [~,hgbase] = EcogExtractHighGammaTrial(sq(ieegBaseCAR(:,iChan,:)), Subject(SN).Experiment.recording.sample_rate, 100, [70 150],baseTimeRange,baseTimeRange,[]);
                
                [m s]=normfit(reshape(log(hgbase(:,51:100)),size(hgbase,1)*length(51:100),1));
                ieegCARHGZ(iChan,:,:)=(log(hgsig)-m)./s;
                ieegBaseCARHGZ(iChan,:,:)=(log(hgbase)-m)./s;
                
                ieegCARHG(iChan,:,:)=hgsig;
                ieegBaseCARHG(iChan,:,:)=hgbase;
                display(iChan);
            end
            
       
            % baseline?
            chanSig={};
            for iChan=1:size(ieegCARHG,1)
                sig1=sq(ieegCARHG(iChan,:,51:size(ieegCARHG,3)-50));
                sig2=repmat(sq(ieegBaseCARHG(iChan,:,51:100)),1,size(sig1,2)./50);
                % sIdx=randperm(300);
                % sig2=sig2(:,sIdx);
                %tic
                [zValsRawAct, pValsRaw, actClust]=timePermCluster(sig1,sig2,1000,1,1.645);
                %toc
                chanSig{iChan}.zValsRawAct=zValsRawAct;
                chanSig{iChan}.pValsRaw=pValsRaw;
                chanSig{iChan}.actClust=actClust;
                chanSig{iChan}.Name=Subject(SN).ChannelInfo(chanIdx(iChan)).Name;
                display(iChan)
            end
            
            if ~exist([DUKEDIR '\Stats\timePerm\'])
                mkdir([DUKEDIR '\Stats\timePerm\'])
            end
            save([DUKEDIR '\Stats\timePerm\' Subject(SN).Name '_' Task.Name '_' ...
                Task.Conds(iC).Name '_' Task.Conds(iC).Field(iF).Name '_' Task.Base.Name '.mat'],'chanSig','ieegCARHG','ieegBaseCARHG','ieegCARHGZ','ieegBaseCARHGZ');
        end
    end
end
