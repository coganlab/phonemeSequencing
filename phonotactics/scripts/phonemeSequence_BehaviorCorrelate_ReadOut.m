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


   
   rMaxOnsetAll=[];
   pMaxOnsetAll=[];
   rMaxDurationAll=[];
   pMaxDurationAll=[];
   rPowerOnsetAll=[];
   pPowerOnsetAll=[];
   rPowerDurationAll=[];
   pPowerDurationAll=[];
   rMaxTimeOnsetAll=[];
   pMaxTimeOnsetAll=[];
   rMaxTimeDurationAll=[];
   pMaxTimeDurationAll=[];

for iSN=1:length(Subject);
    SN=iSN;
   iC=1;
   iF=1;
 load([DUKEDIR '/Stats/BehaviorCorrelates/' Subject(SN).Name '_' ...
     Task.Conds(iC).Name '_' Task.Conds(iC).Field(iF).Name '.mat'])
   rMaxOnsetAll=cat(1,rMaxOnsetAll,rMaxOnset);
   pMaxOnsetAll=cat(1,pMaxOnsetAll,pMaxOnset);
   rMaxDurationAll=cat(1,rMaxDurationAll,rMaxDuration);
   pMaxDurationAll=cat(1,pMaxDurationAll,pMaxDuration);
   rPowerOnsetAll=cat(1,rPowerOnsetAll,rPowerOnset);
   pPowerOnsetAll=cat(1,pPowerOnsetAll,pPowerOnset);
   rPowerDurationAll=cat(1,rPowerDurationAll,rPowerDuration);
   pPowerDurationAll=cat(1,pPowerDurationAll,pPowerDuration);
   rMaxTimeOnsetAll=cat(1,rMaxTimeOnsetAll,rMaxTimeOnset);
   pMaxTimeOnsetAll=cat(1,pMaxTimeOnsetAll,pMaxTimeOnset);
   rMaxTimeDurationAll=cat(1,rMaxTimeDurationAll,rMaxTimeDuration);
   pMaxTimeDurationAll=cat(1,pMaxTimeDurationAll,pMaxTimeDuration);
end

   rMaxOnsetAll(noNameChans)=[];
   pMaxOnsetAll(noNameChans)=[];
   rMaxDurationAll(noNameChans)=[];
   pMaxDurationAll(noNameChans)=[];
   rPowerOnsetAll(noNameChans)=[];
   pPowerOnsetAll(noNameChans)=[];
   rPowerDurationAll(noNameChans)=[];
   pPowerDurationAll(noNameChans)=[];
   rMaxTimeOnsetAll(noNameChans)=[];
   pMaxTimeOnsetAll(noNameChans)=[];
   rMaxTimeDurationAll(noNameChans)=[];
   pMaxTimeDurationAll(noNameChans)=[];
