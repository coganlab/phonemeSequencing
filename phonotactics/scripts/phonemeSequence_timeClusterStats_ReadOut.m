global TASK_DIR
global experiment
global DUKEDIR
global BOX_DIR
global RECONDIR
%BOX_DIR='C:\Users\gcoga\Box';
BOX_DIR='C:\Users\sd355\Box';

%RECONDIR='C:\Users\gcoga\Box\ECoG_Recon';
RECONDIR=[BOX_DIR '\ECoG_Recon'];

Task=[];

Task.Name='Phoneme_Sequencing';
Task.Base.Name='Start';
Task.Base.Epoch='Start';
Task.Base.Time=[-500 0];

Task.Directory=([BOX_DIR '/CoganLab/D_Data/' Task.Name]);
TASK_DIR=Task.Directory;
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
Subject=Subject([1:36]);
Subject(24) = [];
%Subject=Subject([1,3,4,5,6,8:20,24:29]); % precental patients
% badSubjs=[20];
% Subject=Subject(setdiff(1:length(Subject),badSubs));

sigTimeVals={};
subjIdx=1:length(Subject);
timePermDir=[TASK_DIR '\Stats\timePerm\'];
%% Loading significant time points
sigMatAll={};
sigPowerAll={};
for iC=1:length(Task.Conds.Field)
    iC
    for iS=1:length(subjIdx)
        iS
        load([timePermDir Subject(subjIdx(iS)).Name ...
            '_' Subject(subjIdx(iS)).Task '_' Task.Conds.Name ...
            '_' Task.Conds.Field(iC).Name '_' Task.Base.Name '.mat']);
        sigMatTmp=zeros(length(chanSig),length(chanSig{1}.pValsRaw));
        for iChan=1:length(chanSig)
            for iClust=1:length(chanSig{iChan}.actClust.Size)
                if chanSig{iChan}.actClust.Size{iClust}>chanSig{iChan}.actClust.perm99
                    sigMatTmp(iChan,chanSig{iChan}.actClust.Start{iClust}: ...
                        chanSig{iChan}.actClust.Start{iClust}+(chanSig{iChan}.actClust.Size{iClust}-1)) ...
                        = 1;
                end
            end
        end
        if iS==1
            sigMatAll{iC}=sigMatTmp;
           sigPowerAll{iC}=sq(mean(ieegCARHGZ,2));
        else
            sigMatAll{iC}=cat(1,sigMatAll{iC},sigMatTmp);
            sigPowerAll{iC}=cat(1,sigPowerAll{iC},sq(mean(ieegCARHGZ,2)));
        end
    end
end
%% loading all channels
elecNameAllTimePerm = [];
 for iS=1:length(subjIdx)
        iS
        load([timePermDir Subject(subjIdx(iS)).Name ...
            '_' Subject(subjIdx(iS)).Task '_' Task.Conds.Name ...
            '_' Task.Conds.Field(1).Name '_' Task.Base.Name '.mat']);
        elecNameSubj = [];
        for iChan = 1:length(chanSig)
            elecNameSubj{iChan} =  chanSig{iChan}.Name;
        end
        elecNameAllTimePerm= [elecNameAllTimePerm elecNameSubj];
 end
 %% Loading anatomical information
 subj_labels_loc_all = [];
 subj_labels_Name_all = [];
 for iS=1:length(subjIdx)
        iS
    %[subj_labels_Name, subj_labels_loc,subj_labels_WM] = genElecLocationsSubj(Subject(subjIdx(iS)).Name);
    subj_labels_loc_all = [subj_labels_loc_all {Subject(subjIdx(iS)).ChannelInfo.Location}];
    subj_labels_Name_all = [subj_labels_Name_all {Subject(subjIdx(iS)).ChannelInfo.Name}];
 end
 emptyIds = cellfun(@isempty,subj_labels_Name_all);
 subj_labels_Name_all(emptyIds) = [];
 subj_labels_loc_all(emptyIds) = [];
%% 
% There is mismatch between number of channels in the previous and the
% current section. Needs fixing
WMidx=[];%zeros(length(sigMatAll{1}),1);
counter=0;
for iS=1:length(subjIdx);
    WMidx(counter+1:length(Subject(subjIdx(iS)).goodChannels))=1;
    WMvals=intersect(Subject(subjIdx(iS)).WM,Subject(subjIdx(iS)).goodChannels);
    [ii jj]=intersect(Subject(subjIdx(iS)).goodChannels,WMvals);
    WMidx(counter+1:counter+length(Subject(subjIdx(iS)).goodChannels))=0;
    WMidx(counter+jj)=1;
    counter=counter+length(Subject(subjIdx(iS)).goodChannels);
end

    elecNamesAll=[];
    elecLocsAll=[];
for iS=1:length(subjIdx)
    elecNamesTmp=[];
    elecLocsTmp=[];
    for iChan=1:length(Subject(subjIdx(iS)).goodChannels)
        elecNamesTmp{iChan}=Subject(subjIdx(iS)).ChannelInfo(Subject(subjIdx(iS)).goodChannels(iChan)).Name;
        elecLocsTmp{iChan}=Subject(subjIdx(iS)).ChannelInfo(Subject(subjIdx(iS)).goodChannels(iChan)).Location;
    end
    if iS==1
        elecNamesAll=elecNamesTmp;
        elecLocsAll=elecLocsTmp;
    else
        elecNamesAll=cat(2,elecNamesAll,elecNamesTmp);
        elecLocsAll=cat(2,elecLocsAll,elecLocsTmp);
    end
end
%%
% kluge!
noNameChans=[];
counter=0;
for iChan=1:length(elecNameAllTimePerm)
    if isempty(elecNameAllTimePerm{iChan});
        noNameChans(counter+1)=iChan;
        counter=counter+1;
    end
end

% load muscle artifact channels
counter=0;
muscleChannelsAll=[];

for iSN=1:length(Subject);
    load([DUKEDIR '/' Subject(iSN).Name '/muscleChannels.mat'])
    muscleChannels=muscleChannels+counter;
    muscleChannelsAll=cat(2,muscleChannelsAll,muscleChannels);
    counter=counter+length(Subject(iSN).goodChannels);
end
%% Removing null and noisy channels
noNameChans=cat(2,noNameChans,muscleChannelsAll);
sigMatAllClean = sigMatAll;
for iS=1:size(sigMatAllClean,2);
    sigMatAllClean{iS}(noNameChans,:)=[];
    %sigPowerAll{iS}(noNameChans,:)=[];
end
elecNameAllTimePermClean = elecNameAllTimePerm;
elecNameAllTimePermClean(noNameChans)=[];
elecAnatAllTimePermClean = subj_labels_loc_all(ismember(subj_labels_Name_all,elecNameAllTimePermClean));
temporalIds = contains(elecAnatAllTimePermClean,'temporal');
%elecLocsAll(noNameChans)=[];
% WMidx(noNameChans)=[];
% iiNWM=find(WMidx==0);
% iiW=find(WMidx==1);

%% Selecting auditory ids
% audTimeIds = 51:100;
% baseTimeIds = 1:50;
% 
% iiAE=(sum(sigMatAllClean{1}(:,audTimeIds),2)>1);
% iiB = (sum(sigMatAllClean{1}(:,baseTimeIds),2)==0);
% numAudChan(1) = sum(iiAE&iiB);
% figure; imagesc(sigMatAllClean{1}(iiAE&iiB,:));
% for iPer = 2:10
%     iiAE=(sum(sigMatAllClean{1}(:,audTimeIds),2)>length(audTimeIds)*iPer/10);
%     figure; imagesc(sigMatAllClean{1}(iiAE&iiB,:));
%     title([num2str(iPer*100/10) ' % fraction of time activated'])
%     numAudChan(iPer) = sum(iiAE&iiB);
% end

%% Selecting channels that are production, feedback, muscle
% Production baseline = 0, response onset = 1
% Feedback baseline = 0, auditory epoch = 1, pre-response onset = 0 , post-response onset = 1
baseTimeIds = 1:50;
audTimeIds = 51:100;
respTimeIds = 51:150;
preRespTimeIds = 51:100;
postRespTimeIds = 101:150;

baseClustFalse = (sum(sigMatAllClean{1}(:,baseTimeIds),2)==0);
audClustTrue = (sum(sigMatAllClean{1}(:,audTimeIds),2)>1);
preRespClustFalse = (sum(sigMatAllClean{3}(:,preRespTimeIds),2)==0);
postRespClustTrue = (sum(sigMatAllClean{3}(:,postRespTimeIds),2)>1);
respClustTrue = (sum(sigMatAllClean{3}(:,respTimeIds),2)>1);

elecNameProduction = elecNameAllTimePermClean(baseClustFalse...
    &respClustTrue);
elecNameFeedBack = elecNameAllTimePermClean(baseClustFalse&audClustTrue...
    &preRespClustFalse&postRespClustTrue);

feedbackChannelIds = ismember(elecNameProduction,elecNameFeedBack);
muscleChannelIds = ismember(elecNameProduction,elecNameMuscleArtifact);

elecNameProductionClean = elecNameProduction;
elecNameProductionClean(feedbackChannelIds|muscleChannelIds) = [];


%% Selecting channels that are active only during response
% baseline = 0, auditory epoch = 0, pre-response onset = 0 , post-response
% onset = 1
baseTimeIds = 1:50;
audTimeIds = 51:100;
preRespTimeIds = 51:100;
postRespTimeIds = 101:150;
baseClustFalse = (sum(sigMatAllClean{1}(:,baseTimeIds),2)==0);
audClustFalse = (sum(sigMatAllClean{1}(:,audTimeIds),2)==0);
preRespClustFalse = (sum(sigMatAllClean{3}(:,preRespTimeIds),2)==0);
postRespClustTrue = (sum(sigMatAllClean{3}(:,postRespTimeIds),2)>1);
elecNameUtterActive = elecNameAllTimePermClean(baseClustFalse&audClustFalse...
    &preRespClustFalse&postRespClustTrue);

%% Significant activation channel ids for all fields



iiB = (sum(sigMatAllClean{1}(:,baseTimeIds),2)==0);
iiA=(sum(sigMatAllClean{1},2)>1);
iiAE=(sum(sigMatAllClean{1}(:,51:100),2)>1);
iiAE2=(sum(sigMatAllClean{1}(:,51:75),2)>1);

iiG=(sum(sigMatAllClean{2},2)>1);
iiR=(sum(sigMatAllClean{3}(:,51:150),2)>1);

%iiR2=find(sum(sigMatAll{3}(:,75:125),2)>1);
iiS=(sum(sigMatAllClean{4},2)>1);
iiD=(sum(sigMatAllClean{1}(:,150:200),2)>1);


iiSM=intersect(iiA,iiR);
iiA2=setdiff(iiA,iiSM);
iiR2=setdiff(iiR,iiSM);





iiRD=intersect(iiR,iiD);
%%
iiDNW=intersect(iiD,iiNWM);
iiRNW=intersect(iiR,iiNWM);
iiSMNW=intersect(iiSM,iiNWM);
iiA2NW=intersect(iiA2,iiNWM);
iiR2NW=intersect(iiR2,iiNWM);
iiRDNW=intersect(iiRNW,iiDNW);

for iChan=1:size(sigMatAllClean{3},1);
    [ii jj]=find(sigMatAllClean{3}(iChan,:)==1);
    if ~isempty(jj)
        iiOnset(iChan)=jj(1);
    else
        iiOnset(iChan)=NaN;
    end
end
    
iiEarly=find(iiOnset<90);
iiLate=find(iiOnset>=90);
%%
grouping_idx=zeros(length(elecNamesAll),1);
 % grouping_idx(iiA2NW)=2;% green
  grouping_idx(iiRD)=1; % red
 %   grouping_idx(intersect(iiR2NW,iiEarly))=2; % green

%  grouping_idx(setdiff(iiSMNW,iiDNW))=3; % yellow
% grouping_idx(intersect(iiR2NW,iiEarly))=4; % blue
% %grouping_idx(intersect(iiR2NW,iiLate))=2;

%grouping_idx(iiRDNWEarly)=1;
%grouping_idx(iiRNDNWEarly)=4;
%grouping_idx(labelChanBroca)=1;
cfg=[];
cfg.hemisphere='l';
plot_subjs_on_average_grouping(elecNamesAll, grouping_idx, 'fsaverage', cfg)


% labels?
labelChan=[];
labelChanName={};
%label='Plan_temp';
%label='precentral';
%label='postcentral';
%label='superiortemporal'
%label='supram'
%label='Triang';
label='Operc';
%label='S_central';
counter=0;
for iChan=1:length(elecLocsAll);
    if contains(elecLocsAll{iChan},label,'IgnoreCase',1);
        labelChan(counter+1)=iChan;
        labelChanName{counter+1}=elecNamesAll{iChan};
        counter=counter+1;
    end
end

DNWLocs=elecLocsAll(iiRNDNWEarly);
[uni,~,idx]=unique(DNWLocs);
figure;
hist(idx,unique(idx));