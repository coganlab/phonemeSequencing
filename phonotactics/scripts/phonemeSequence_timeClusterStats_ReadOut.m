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
% Subject=Subject(1:36);
% Subject(24) = [];
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
%            sigPowerAll{iC}=sq(mean(ieegCARHGZ,2));
        else
            sigMatAll{iC}=cat(1,sigMatAll{iC},sigMatTmp);
%             sigPowerAll{iC}=cat(1,sigPowerAll{iC},sq(mean(ieegCARHGZ,2)));
        end
    end
end

% loading all channels
elecNameAllTimePerm = [];
 for iS=1:length(subjIdx)
        iS
        load([timePermDir Subject(subjIdx(iS)).Name ...
            '_' Subject(subjIdx(iS)).Task '_' Task.Conds.Name ...
            '_' Task.Conds.Field(1).Name '_' Task.Base.Name '.mat']);
        elecNameSubj = [];
        assert(length(chanSig)==length(channelNames),'Channel dimension mismatch');
%         for iChan = 1:length(chanSig)
%             elecNameSubj{iChan} =  channelNames;
%         end
        elecNameAllTimePerm= [elecNameAllTimePerm channelNames];
 end

 emptyIds = cellfun(@isempty,elecNameAllTimePerm);
 elecNameAllTimePerm(emptyIds) = [];
 for iS=1:size(sigMatAll,2)
    sigMatAll{iS}(emptyIds,:)=[];
    %sigPowerAll{iS}(noNameChans,:)=[];
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

%% Loading Muscle Channels
muscleChannelsAll = []
for iSN=1:length(Subject);
    Subject(iSN).Name
    if exist([DUKEDIR '/' Subject(iSN).Name '/muscleChannelWavelet.mat'])
     load([DUKEDIR '/' Subject(iSN).Name '/muscleChannelWavelet.mat'])
     try
        muscleChannelsAll=cat(2,muscleChannelsAll,(muscleChannel));
     catch
        muscleChannelsAll=cat(2,muscleChannelsAll,(muscleChannel)');
     end
    end
 %   muscleChannels=muscleChannels+counter;
%     muscleChannelsAll=cat(2,muscleChannelsAll,muscleChannels);
%     counter=counter+length(Subject(iSN).goodChannels);
end
muscleChannelsAll =  muscleChannelsAll(~cellfun('isempty',muscleChannelsAll));
%% Removing Unknown & muscle channels
% 

muscleChannelIds = ismember(elecNameAllTimePerm,muscleChannelsAll);
unknownChannels = subj_labels_Name_all(contains(subj_labels_loc_all,'Unknown'));
unknownChannelIds = ismember(elecNameAllTimePerm,unknownChannels);
elecNameAllTimePermClean =  elecNameAllTimePerm(~(muscleChannelIds|unknownChannelIds));
for iS=1:size(sigMatAll,2)
    iS
    sigMatClean{iS}=sigMatAll{iS}(~(muscleChannelIds|unknownChannelIds),:);
    %sigPowerAll{iS}(noNameChans,:)=[];
end


%% Selecting channels that are feedback
% baseline = 0, auditory epoch = 1, pre-response onset = 0 , post-response
% onset = 1
timeAuditory = linspace(-0.5,1.5,200);
timeResponse = linspace(-1,1,200);
baseTimeIds = 1:50;
audTimeIds = timeAuditory>0&timeAuditory<0.2;
preRespTimeIds = timeResponse>-0.5&timeResponse<-0.1;
postRespTimeIds = timeResponse>=-0.1;
baseClustFalse = (sum(sigMatClean{4}(:,baseTimeIds),2)==0);
audClustTrue = (sum(sigMatClean{1}(:,audTimeIds),2)>1);
preRespClustFalse = (sum(sigMatClean{3}(:,preRespTimeIds),2)==0);
postRespClustTrue = (sum(sigMatClean{3}(:,postRespTimeIds),2)>1);
feedbackCond = baseClustFalse&audClustTrue...
    &preRespClustFalse&postRespClustTrue;
elecNameFeedBack_100 = elecNameAllTimePermClean(feedbackCond);
elecNameFeedBack_100 =  elecNameFeedBack_100(~cellfun('isempty',elecNameFeedBack_100));

elecAnatFeedBack_100 = subj_labels_loc_all(ismember(subj_labels_Name_all,elecNameFeedBack_100));
anat2select = ["temporal","White","sts"];
elecNameFeedBack_100_temporal = elecNameFeedBack_100(contains(elecAnatFeedBack_100,anat2select));


for iS=1:size(sigMatClean,2)
    sigMatFeedbackTemporal{iS}=sigMatClean{iS}(ismember(elecNameAllTimePermClean,elecNameFeedBack_100),:);
    %sigPowerAll{iS}(noNameChans,:)=[];
end

figure;
subplot(1,2,1)
imagesc(timeAuditory,[],sigMatFeedbackTemporal{1});
xlabel('Time from Auditory onset (s)');
subplot(1,2,2)
imagesc(timeResponse,[],sigMatFeedbackTemporal{3});
xlabel('Time from Response onset (s)');
sgtitle('Feedback channels')




%% Selecting channels that are production, feedback, muscle
% Production baseline = 0, response onset = 1
% Feedback baseline = 0, auditory epoch = 1, pre-response onset = 0 , post-response onset = 1
baseTimeIds = 1:50;
audTimeIds = find(timeAuditory>0&timeAuditory<0.2);
respTimeIds = timeResponse>-0.5&timeResponse<0.5;
preRespTimeIds = find(timeResponse>-0.5&timeResponse<-0.1);
postRespTimeIds = find(timeResponse>-0.1&timeResponse<0.5);

baseClustFalse = (sum(sigMatClean{4}(:,baseTimeIds),2)==0);
audBaseClustFalse = (sum(sigMatClean{1}(:,baseTimeIds),2)==0);
audClustTrue = (sum(sigMatClean{1}(:,audTimeIds),2)>1);
audClustFalse = (sum(sigMatClean{1}(:,audTimeIds),2)==0);
% preRespClustFalse = (sum(sigMatAll{3}(:,preRespTimeIds),2)==0);
% postRespClustTrue = (sum(sigMatAll{3}(:,postRespTimeIds),2)>1);
respClustTrue = (sum(sigMatClean{3}(:,respTimeIds),2)>1);

elecNameAudTrueProdTrue = elecNameAllTimePermClean(baseClustFalse&...
    audClustTrue&respClustTrue);
elecNameAudTrueProdTrue =  elecNameAudTrueProdTrue(~cellfun('isempty',elecNameAudTrueProdTrue));
%elecNameAudTrueProdTrueExcludeMuscle = elecNameAudTrueProdTrue(~ismember(elecNameAudTrueProdTrue,muscleChannelsAll));
elecNameAudTrueProdTrueExcludeFeedback = elecNameAudTrueProdTrue(~ismember(elecNameAudTrueProdTrue,elecNameFeedBack_100));

tempLobeAnat2select = ["parahippocampal","fusiform","entorhinal","Cerebellum","lingual","Ventricle"];
elecAnatAudTrueProdTrueExcludeFeedback = subj_labels_loc_all(ismember(subj_labels_Name_all,elecNameAudTrueProdTrueExcludeFeedback));
temporalLobeChannels = contains(elecAnatAudTrueProdTrueExcludeFeedback,tempLobeAnat2select);
%elecNameAudTrueProdTrueNoTemporalLobe = elecNameAudTrueProdTrueExcludeFeedback(~temporalLobeChannels);



for iS=1:size(sigMatClean,2)
    sigMatAudTrueProdTrueNoFeedback{iS}=sigMatClean{iS}(ismember(elecNameAllTimePermClean,elecNameAudTrueProdTrueExcludeFeedback),:);
    %sigPowerAll{iS}(noNameChans,:)=[];
end

figure;
subplot(1,2,1)
imagesc(timeAuditory,[],sigMatAudTrueProdTrueNoFeedback{1});
xlabel('Time from Auditory onset (s)');
subplot(1,2,2)
imagesc(timeResponse,[],sigMatAudTrueProdTrueNoFeedback{3});
xlabel('Time from Response onset (s)');
sgtitle('Auditory True + Production true (No feedback)')


% elecAnatAudTrueProdTrueExcludeMuscleExcludeFeedback = subj_labels_loc_all(ismember(subj_labels_Name_all,elecNameAudTrueProdTrueExcludeMuscleExcludeFeedback));
% elecNameAudTrueProdTrueTemporal = elecNameAudTrueProdTrueExcludeMuscle(contains(elecAnatAudTrueProdTrueExcludeMuscleExcludeFeedback,'temporal'));
% 
% for iS=1:size(sigMatAll,2)
%     sigMatFeedbackAudTrueProdTrueTemporal{iS}=sigMatAll{iS}(ismember(elecNameAllTimePerm,elecNameAudTrueProdTrueTemporal),:);
%     %sigPowerAll{iS}(noNameChans,:)=[];
% end
% 
% figure;
% subplot(1,2,1)
% imagesc(timeAuditory,[],sigMatFeedbackAudTrueProdTrueTemporal{1});
% xlabel('Time from Auditory onset (s)');
% subplot(1,2,2)
% imagesc(timeResponse,[],sigMatFeedbackAudTrueProdTrueTemporal{3});
% xlabel('Time from Response onset (s)');
% sgtitle('Auditory True + Production true (Temporal)')
% 
% preRespTimeFeedbackIds = find(timeResponse<=-0.1);
% postRespTimeFeedbackIds = find(timeResponse>-0.1);
% elecNameAudTrueProdTrueTemporalFeedback = elecNameAudTrueProdTrueTemporal(sum(sigMatFeedbackAudTrueProdTrueTemporal{3}(:,preRespTimeIds),2)==0&sum(sigMatFeedbackAudTrueProdTrueTemporal{3}(:,postRespTimeIds),2)>1)';
% 


% sigDiff = (diff(sigMatFeedbackAudTrueProdTrueTemporal{3}'));
% for iChan = 1:size(sigDiff,2)
%     respOnsId = find(sigDiff(:,iChan)==1,1);
%     if(isempty(respOnsId))
%         respOnsId = 1;
%     end
%     respOnsAudTrueTemporal(iChan) = timeResponse(respOnsId);
% end



%%


% elecAnatAudTrueProdTrue = subj_labels_loc_all(ismember(subj_labels_Name_all,elecNameAudTrueProdTrue));
% 
% elecNameAudTrueProdTrueTemporal = elecNameAudTrueProdTrue(contains(elecAnatAudTrueProdTrue,'temporal'));
% elecNameAudTrueProdTrueTemporalExcludeFeedback = elecNameAudTrueProdTrueTemporal(~ismember(elecNameAudTrueProdTrueTemporal,elecNameFeedBack_100));


elecNameAudFalseProdTrue = elecNameAllTimePermClean(baseClustFalse&...
    audClustFalse&respClustTrue);
elecNameAudFalseProdTrue =  elecNameAudFalseProdTrue(~cellfun('isempty',elecNameAudFalseProdTrue));
%elecNameAudFalseProdTrueExcludeMuscle = elecNameAudFalseProdTrue(~ismember(elecNameAudFalseProdTrue,muscleChannelsAll));

for iS=1:size(sigMatClean,2)
    sigMatFeedbackAudFalseProdTrue{iS}=sigMatClean{iS}(ismember(elecNameAllTimePermClean,elecNameAudFalseProdTrue),:);
    %sigPowerAll{iS}(noNameChans,:)=[];
end

figure;
subplot(1,2,1)
imagesc(timeAuditory,[],sigMatFeedbackAudFalseProdTrue{1});
xlabel('Time from Auditory onset (s)');
subplot(1,2,2)
imagesc(timeResponse,[],sigMatFeedbackAudFalseProdTrue{3});
xlabel('Time from Response onset (s)');
sgtitle('Auditory False + Production true ')



elecAnatAudFalseProdTrue = subj_labels_loc_all(ismember(subj_labels_Name_all,elecNameAudFalseProdTrue));

tempAnat2select = ["temporal","sts"];
temporalChannels = contains(elecAnatAudFalseProdTrue,tempAnat2select);
tempLobeAnat2select = ["Hippocampus","Amygdala","parahippocampal","fusiform","entorhinal","Cerebellum","lingual","Ventricle"];
temporalLobeChannels = contains(elecAnatAudFalseProdTrue,tempLobeAnat2select);

elecNameAudFalseProdTrueTemporal = elecNameAudFalseProdTrue(temporalChannels|temporalLobeChannels);

for iS=1:size(sigMatClean,2)
    sigMatFeedbackAudFalseProdTrueTemporal{iS}=sigMatClean{iS}(ismember(elecNameAllTimePermClean,elecNameAudFalseProdTrueTemporal),:);
    %sigPowerAll{iS}(noNameChans,:)=[];
end

figure;
subplot(1,2,1)
imagesc(timeAuditory,[],sigMatFeedbackAudFalseProdTrueTemporal{1});
xlabel('Time from Auditory onset (s)');
subplot(1,2,2)
imagesc(timeResponse,[],sigMatFeedbackAudFalseProdTrueTemporal{3});
xlabel('Time from Response onset (s)');
sgtitle('Auditory False + Production true (only temporal)')


sigDiff = (diff(sigMatFeedbackAudFalseProdTrueTemporal{3}'));
for iChan = 1:size(sigDiff,2)
    respOnsId = find(sigDiff(:,iChan)==1,1);
    respOns(iChan) = timeResponse(respOnsId);
end


elecNameAudFalseProdTrueNoTemporal = elecNameAudFalseProdTrue(~(temporalChannels|temporalLobeChannels));

for iS=1:size(sigMatClean,2)
    sigMatFeedbackAudFalseProdTrueNoTemporal{iS}=sigMatClean{iS}(ismember(elecNameAllTimePermClean,elecNameAudFalseProdTrueNoTemporal),:);
    %sigPowerAll{iS}(noNameChans,:)=[];
end

figure;
subplot(1,2,1)
imagesc(timeAuditory,[],sigMatFeedbackAudFalseProdTrueNoTemporal{1});
xlabel('Time from Auditory onset (s)');
subplot(1,2,2)
imagesc(timeResponse,[],sigMatFeedbackAudFalseProdTrueNoTemporal{3});
xlabel('Time from Response onset (s)');
sgtitle('Auditory False + Production true (Temporal Excluded) ')

motorChannel = elecNameAudFalseProdTrueNoTemporal;

motorAnat = subj_labels_loc_all(ismember(subj_labels_Name_all,motorChannel));
motorChannelGM = motorChannel(~contains(motorAnat,'White'));

%elecNameAudFalseProdTrueTemporalExcludeMuscle = elecNameAudFalseProdTrueTemporal(~ismember(elecNameAudFalseProdTrueTemporal,muscleChannelsAll));

% elecNameAudFalseProdTrueFrontal = elecNameAudFalseProdTrue(contains(elecAnatAudFalseProdTrue,'frontal'));
% elecNameAudFalseProdTrueFrontalExcludeMuscle = elecNameAudFalseProdTrueFrontal(~ismember(elecNameAudFalseProdTrueFrontal,muscleChannelsAll));
% elecNameFeedBack = elecNameAllTimePerm(baseClustFalse&audClustTrue...
%     &preRespClustFalse&postRespClustTrue);
% elecNameFeedBack =  elecNameFeedBack(~cellfun('isempty',elecNameFeedBack));

% feedbackChannelIds = ismember(elecNameProduction,elecNameFeedBack_100);
% muscleChannelIds = ismember(elecNameProduction,muscleChannelsAll);
% elecNameProductionClean = elecNameProduction;
% elecNameProductionClean(feedbackChannelIds|muscleChannelIds) = [];


%% Selecting channels that are production, feedback, muscle
% Production baseline = 0, response onset = 1
% Feedback baseline = 0, auditory epoch = 1, pre-response onset = 0 , post-response onset = 1
baseTimeIds = 1:50;

respTimeIds = find(timeResponse>-0.5&timeResponse<0.5);


baseClustFalse = (sum(sigMatAll{4}(:,baseTimeIds),2)==0);

% preRespClustFalse = (sum(sigMatAll{3}(:,preRespTimeIds),2)==0);
% postRespClustTrue = (sum(sigMatAll{3}(:,postRespTimeIds),2)>1);
respClustTrue = (sum(sigMatAll{3}(:,respTimeIds),2)>1);

elecNameProduction = elecNameAllTimePerm(baseClustFalse...
    &respClustTrue);
elecNameProduction =  elecNameProduction(~cellfun('isempty',elecNameProduction));
% elecNameFeedBack = elecNameAllTimePerm(baseClustFalse&audClustTrue...
%     &preRespClustFalse&postRespClustTrue);
% elecNameFeedBack =  elecNameFeedBack(~cellfun('isempty',elecNameFeedBack));

feedbackChannelIds = ismember(elecNameProduction,elecNameFeedBack_100_temporal);
muscleChannelIds = ismember(elecNameProduction,muscleChannelsAll);
temporalChannelAudTrueProdTrueIds = ismember(elecNameProduction,elecNameAudTrueProdTrueTemporal);
temporalChannelAudFalseProdTrueIds = ismember(elecNameProduction,elecNameAudFalseProdTrueTemporal);
elecNameProductionClean = elecNameProduction;
elecNameProductionClean(feedbackChannelIds|muscleChannelIds|temporalChannelAudTrueProdTrueIds|temporalChannelAudFalseProdTrueIds) = [];


% elecNameProduction = elecNameAllTimePerm(baseClustFalse&audBaseClustFalse...
%     &respClustTrue);
% elecNameProduction =  elecNameProduction(~cellfun('isempty',elecNameProduction));
% elecNameFeedBack = elecNameAllTimePerm(baseClustFalse&audBaseClustFalse&audClustTrue...
%     &preRespClustFalse&postRespClustTrue);
% elecNameFeedBack =  elecNameFeedBack(~cellfun('isempty',elecNameFeedBack));
% 
% feedbackChannelIds = ismember(elecNameProduction,elecNameFeedBack);
% elecNameProductionCleanAudBase = elecNameProduction;
% elecNameProductionCleanAudBase(feedbackChannelIds) = [];


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