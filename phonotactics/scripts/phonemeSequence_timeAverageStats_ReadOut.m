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


Task.Directory=([BOX_DIR '/CoganLab/D_Data/' Task.Name]);
TASK_DIR=Task.Directory;
DUKEDIR=TASK_DIR;

Task.Conds(1).Name='All';
Task.Conds(1).Field(1).Name='AuditorywDelay';
Task.Conds(1).Field(1).Epoch='Auditory';

Task.Conds(1).Field(2).Name='DelaywGo';
Task.Conds(1).Field(2).Epoch='Go';

Task.Conds(1).Field(3).Name='Response';
Task.Conds(1).Field(3).Epoch='ResponseStart';



Subject = popTaskSubjectData(Task);
% Subject=Subject(1:36);
% Subject(24) = [];
%Subject=Subject([1,3,4,5,6,8:20,24:29]); % precental patients
% badSubjs=[20];
% Subject=Subject(setdiff(1:length(Subject),badSubs));

sigTimeVals={};
subjNames = extractfield(Subject,'Name');
% subjList = {'D76','D80','D81','D84','D85','D86','D88','D92','D93','D94'...
%      'D95','D96','D103'};
subjList = extractSubjectPerChannel(elecNameProductionClean);
subjIdx = find(ismember(subjNames,subjList));

timeAveDir=[TASK_DIR '\Stats\timeGroup\'];
%% Loading significant time points
sigMatAll={};
sigPowerAll={};
for iC=1:3
    iC
    for iS=1:length(subjIdx)
        iS
        load([timeAveDir Subject(subjIdx(iS)).Name ...
            '_' Subject(subjIdx(iS)).Task '_' Task.Conds.Name ...
            '_' Task.Conds.Field(iC).Name '_' Task.Base.Name '.mat']);
        sigMatTmp=zeros(1,length(pChan));
%         for iChan=1:length(chanSig)
%             for iClust=1:length(chanSig{iChan}.actClust.Size)
%                 if chanSig{iChan}.actClust.Size{iClust}>chanSig{iChan}.actClust.perm99
%                     sigMatTmp(iChan,chanSig{iChan}.actClust.Start{iClust}: ...
%                         chanSig{iChan}.actClust.Start{iClust}+(chanSig{iChan}.actClust.Size{iClust}-1)) ...
%                         = 1;
%                 end
%             end
%         end
        [~,p_masked] = fdr(pChan,0.05);
        sigMatTmp = p_masked;
        if iS==1
            sigMatAll{iC}=sigMatTmp;
%            sigPowerAll{iC}=sq(mean(ieegCARHGZ,2));
        else
            sigMatAll{iC}=cat(2,sigMatAll{iC},sigMatTmp);
%             sigPowerAll{iC}=cat(1,sigPowerAll{iC},sq(mean(ieegCARHGZ,2)));
        end
    end
end

% loading all channels
elecNameAllTimeAveNew = [];
 for iS=1:length(subjIdx)
        iS
        load([timeAveDir Subject(subjIdx(iS)).Name ...
            '_' Subject(subjIdx(iS)).Task '_' Task.Conds.Name ...
            '_' Task.Conds.Field(1).Name '_' Task.Base.Name '.mat']);
        elecNameSubj = [];
        assert(length(pChan)==length(channelNames),'Channel dimension mismatch');
%         for iChan = 1:length(chanSig)
%             elecNameSubj{iChan} =  channelNames;
%         end
        elecNameAllTimeAveNew= [elecNameAllTimeAveNew channelNames];
 end

 emptyIds = cellfun(@isempty,elecNameAllTimeAveNew);
 elecNameAllTimeAveNew(emptyIds) = [];
 for iS=1:size(sigMatAll,2)
    sigMatAll{iS}(emptyIds)=[];
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
        muscleChannelsAll=cat(1,muscleChannelsAll,(muscleChannel));
     catch
        muscleChannelsAll=cat(1,muscleChannelsAll,(muscleChannel)');
     end
    end
 %   muscleChannels=muscleChannels+counter;
%     muscleChannelsAll=cat(2,muscleChannelsAll,muscleChannels);
%     counter=counter+length(Subject(iSN).goodChannels);
end
muscleChannelsAll =  muscleChannelsAll(~cellfun('isempty',muscleChannelsAll));
%% Removing Unknown & muscle channels
% 

muscleChannelIds = ismember(elecNameAllTimeAveNew,muscleChannelsAll);
unknownChannels = subj_labels_Name_all(contains(subj_labels_loc_all,'Unknown'));
unknownChannelIds = ismember(elecNameAllTimeAveNew,unknownChannels);
elecNameAllTimePermCleanNew =  elecNameAllTimeAveNew(~(muscleChannelIds|unknownChannelIds));
for iS=1:size(sigMatAll,2)
    iS
    sigMatClean{iS}=sigMatAll{iS}(~(muscleChannelIds|unknownChannelIds));
    %sigPowerAll{iS}(noNameChans,:)=[];
end


%% Selecting production channels

elecNameProductionNew = elecNameAllTimePermCleanNew(sigMatClean{3})';
elecNameDelay = elecNameAllTimePermCleanNew(sigMatClean{2})';
elecNameAuditory = elecNameAllTimePermCleanNew(sigMatClean{1})';
elecNameProductionInfo = extractChannelLocation(Subject,elecNameProductionNew);
elecNameProductionGM = elecNameProduction(~contains({elecNameProductionInfo.Location},'White'));
elecNameProductionSM = elecNameProduction(contains({elecNameProductionInfo.Location},["central"]));
elecNameProductionTemp = elecNameProduction(contains({elecNameProductionInfo.Location},["temporal","sts"]));
elecNameProductionIfg = elecNameProduction(contains({elecNameProductionInfo.Location},["triang","opercul"]));
elecNameProductionFront = elecNameProduction(contains({elecNameProductionInfo.Location},["middlefrontal"]));
elecNameProductionPariet = elecNameProduction(contains({elecNameProductionInfo.Location},["supramarginal","inferiorparietal"]));


%% Loading Normalized High-Gamma for an ROI and corresponding trialInfo

 fieldEpoch = 'ResponseStart';


selectRoi ={};



respTimeThresh = 0.1;
respDurThresh = 1.5;
timeEpoch = [-2  2];


subSelectElecs = elecNameProductionNew;

ieegHGStruct = extractHGDataWithROI(Subject,baseName = 'Start',...
    Epoch = fieldEpoch, roi = selectRoi, Time=timeEpoch,respTimeThresh=respTimeThresh,...
    subsetElec=subSelectElecs,remWMchannels=false,normType=1,fDown = 100,respDurThresh=respDurThresh);

% Remove empty subjects
emptyIds = [];

trialInfoStruct = [];
for iSubject = 1:length(Subject)
    
    if(isempty(ieegHGStruct(iSubject).ieegHGNorm))
        trialInfoStruct(iSubject).phonemeTrial = [];
        emptyIds = [emptyIds iSubject];
    else
      trialInfoStruct(iSubject).phonemeTrial = phonemeSequenceTrialParser(ieegHGStruct(iSubject).trialInfo);   
    end

end

ieegHGStruct(emptyIds) = [];
trialInfoStruct(emptyIds) = [];

maxTrial = [];
for iSubj = 1:length(ieegHGStruct)
    chanSubj(iSubj) = length(ieegHGStruct(iSubj).channelName);
    trialSubj(iSubj) = length(ieegHGStruct(iSubj).trialInfo);
    if(trialSubj(iSubj)>=156)
        maxTrial = [maxTrial iSubj];
    end
end
%%

[ieegStructPooled,phonemeTrialPooled,channelNamePooled] = poolChannelWithMaxTrial(ieegHGStruct,trialInfoStruct);

ieegMean = squeeze(nanmean(ieegStructPooled.data,2));
absDiffIeeg = (diff(ieegMean'));

%chanNoiseId = (min(absDiffIeeg(105:125,:))>0.25)|(max(absDiffIeeg(368:382,:))>1)|(max(absDiffIeeg(1:399,:))>0.25)
chanNoiseId = (min(absDiffIeeg(1:399,:))<-0.25);
chanNoise = channelNamePooled(chanNoiseId);
elecNameProductionNewClean = elecNameProductionNew(~ismember(elecNameProductionNew,chanNoise));
