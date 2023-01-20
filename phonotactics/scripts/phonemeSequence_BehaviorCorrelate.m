global TASK_DIR
global experiment
global DUKEDIR
global BOX_DIR
%BOX_DIR='C:\Users\gcoga\Box';
BOX_DIR='H:\Box Sync';
%RECONDIR='C:\Users\gcoga\Box\ECoG_Recon';
RECONDIR=[BOX_DIR '\ECoG_Recon'];

Task=[];

Task.Name='Phoneme_Sequencing';


TASK_DIR=([BOX_DIR '/CoganLab/D_Data/' Task.Name]);
DUKEDIR=TASK_DIR;

Task.Conds(1).Name='All';
% Task.Conds(1).Field(1).Name='AuditorywDelay';
% Task.Conds(1).Field(1).Epoch='Auditory';
% Task.Conds(1).Field(1).Time=[-500 1500];
Task.Conds(1).Field(1).Name='DelaywGo';
Task.Conds(1).Field(1).Epoch='Go';
Task.Conds(1).Field(1).Time=[-1000 1500];
% Task.Conds(1).Field(3).Name='Response';
% Task.Conds(1).Field(3).Epoch='ResponseStart';
% Task.Conds(1).Field(3).Time=[-1000 1000];
% Task.Conds(1).Field(4).Name='StartCue';
% Task.Conds(1).Field(4).Epoch='Start';
% Task.Conds(1).Field(4).Time=[-500 1000];
Task.Base.Name='Start';
Task.Base.Epoch='Start';
Task.Base.Time=[-500 0];

Subject = popTaskSubjectData(Task.Name);
Subject=Subject(1:27);

for iSN=1:length(Subject);
    SN=iSN;
       ieegB=trialIEEG(Subject(SN).Trials,Subject(SN).goodChannels,...
  Task.Base.Epoch,[Task.Base.Time(1)-500 Task.Base.Time(2)+500]);
ieegBCAR=ieegB-mean(ieegB,2);
%     for iC=1:length(Task.Conds)
%         for iF=1:length(Task.Conds.Field)
            iC=1;iF=1;
    ieeg=trialIEEG(Subject(SN).Trials,Subject(SN).goodChannels,...
    Task.Conds(iC).Field(iF).Epoch,...
    [Task.Conds(iC).Field(iF).Time(1)-500 Task.Conds(iC).Field(iF).Time(2)+500]);    
   ieegCAR=ieeg-mean(ieeg,2);
   
   ResponseOnset=[];
   ResponseDuration=[];
   for iTrials=1:length(Subject(SN).Trials)
       if ~isempty(Subject(SN).Trials(iTrials).ResponseStart)
           ResponseOnset(iTrials)=(Subject(SN).Trials(iTrials).ResponseStart-Subject(SN).Trials(iTrials).Go)./30000;
           ResponseDuration(iTrials)=(Subject(SN).Trials(iTrials).ResponseEnd-Subject(SN).Trials(iTrials).ResponseStart)./30000;
       else
           ResponseOnset(iTrials)=NaN;
           ResponseDuration(iTrials)=NaN;
       end
   end
   
   rMaxOnset=zeros(size(ieegCAR,2),1);
   pMaxOnset=zeros(size(ieegCAR,2),1);
   rMaxDuration=zeros(size(ieegCAR,2),1);
   pMaxDuration=zeros(size(ieegCAR,2),1);
   rPowerOnset=zeros(size(ieegCAR,2),1);
   pPowerOnset=zeros(size(ieegCAR,2),1);
   rPowerDuration=zeros(size(ieegCAR,2),1);
   pPowerDuration=zeros(size(ieegCAR,2),1);
   rMaxTimeOnset=zeros(size(ieegCAR,2),1);
   pMaxTimeOnset=zeros(size(ieegCAR,2),1);
   rMaxTimeDuration=zeros(size(ieegCAR,2),1);
   pMaxTimeDuration=zeros(size(ieegCAR,2),1);
   
for iChan=1:size(ieegCAR,2)
     [bh,bh2]=EcogExtractHighGammaTrial(sq(ieegBCAR(:,iChan,:)),...
        Subject(SN).Experiment.recording.sample_rate,200,[70 150],...
        [(Task.Base.Time(1)-500)/1000 (Task.Base.Time(2)+500)/1000],...
        [(Task.Base.Time(1))/1000 (Task.Base.Time(2))/1000],[],[]);
[m s]=normfit(bh2(:)); 
    
    [dh,dh2]=EcogExtractHighGammaTrial(sq(ieegCAR(:,iChan,:)),...
        Subject(SN).Experiment.recording.sample_rate,200,[70 150],...
        [(Task.Conds(iC).Field(iF).Time(1)-500)/1000 (Task.Conds(iC).Field(iF).Time(2)+500)/1000],...
        [(Task.Conds(iC).Field(iF).Time(1))/1000 (Task.Conds(iC).Field(iF).Time(2))/1000],[m s],1);
    dh2=dh2(:,101:size(dh2,2)-100);
    
    iiNaN=isnan(ResponseOnset);
    iiGood=find(iiNaN==0);
    
 
    % Max
    [Y I]=max(dh2(iiGood,:),[],2);
    
    [r p]=corrcoef(ResponseOnset(iiGood),max(dh2(iiGood,:),[],2));
    rMaxOnset(iChan)=r(1,2);
    pMaxOnset(iChan)=p(1,2);

    [r p]=corrcoef(ResponseDuration(iiGood),max(dh2(iiGood,:),[],2));
    rMaxDuration(iChan)=r(1,2);
    pMaxDuration(iChan)=p(1,2); 
    
    [r p]=corrcoef(ResponseOnset(iiGood),I);
    rMaxTimeOnset(iChan)=r(1,2);
    pMaxTimeOnset(iChan)=p(1,2);
    
   [r p]=corrcoef(ResponseDuration(iiGood),I);
    rMaxTimeDuration(iChan)=r(1,2);
    pMaxTimeDuration(iChan)=p(1,2);
    % mean power
   [r p]=corrcoef(ResponseOnset(iiGood),mean(dh2(iiGood,:).^2,2));
    rPowerOnset(iChan)=r(1,2);
    pPowerOnset(iChan)=p(1,2);
    [r p]=corrcoef(ResponseDuration(iiGood),mean(dh2(iiGood,:).^2,2));
    rPowerDuration(iChan)=r(1,2);
    pPowerDuration(iChan)=p(1,2); 
    
end
 save([DUKEDIR '/Stats/BehaviorCorrelates/' Subject(SN).Name '_' ...
     Task.Conds(iC).Name '_' Task.Conds(iC).Field(iF).Name '.mat']...
     , 'rMax*', 'rPower*','pMax*','pPower*');
end