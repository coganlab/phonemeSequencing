addpath(genpath('E:\Box Sync\ECoG_Recon\matlab_code\'));
global DUKEDIR
DUKEDIR = 'E:\Box Sync\CoganLab\D_Data\Phoneme_Sequencing\';
dLabels = dir(DUKEDIR);
dLabels = dLabels(3:end);

tw = [-3 2]; % time window
etw = [-2.5 1.5];
etwG = [-1 1];% epoch time window
prtw = [-2.5 -2]; % preonset time window
pstw = [-0.25 0.25]; % postonset time window
gammaF = [70 150]; % frequency in Hz
%%
for iSubject = 1:10
    
    d = []; ieegCarResp = []; ieegCarImpFilt = []; ieegGamma = []; ieegSplit = [];
    Experiment = loadExperiment(dLabels(iSubject).name);
    fsD = Experiment.recording.sample_rate;
    Trials = dbTrials(dLabels(iSubject).name,Experiment.recording.recording_day,'Speech_OvertMimeMove');
    allChannels = string({Experiment.channels.name});
    trialFiles = strcat('\',Experiment.recording.recording_day,'\mat\trialInfo.mat');
    load([DUKEDIR '\' dLabels(iSubject).name '\' trialFiles])
    %% Recon Visualization
%     iSubject = 1;
%     dLabels(iSubject).name = 'D2';
%     cfg = [];
%     % cfg.alpha = 0.4;
%     % cfg.use_brainshifted = 1;
%     handles = plot_subj(dLabels(iSubject).name, cfg);
    
channelName = {Subject(iSubject).ChannelInfo.Location};
motorChan = contains(channelName,'precentral');
sensoryChan = contains(channelName,'postcentral');
ifgChan = contains(channelName,'opercularis');

anatChan = motorChan|sensoryChan|ifgChan;

%%
    inUnitFile = allChannels';
    gridIndex = zeros(size(inUnitFile));
    ieegCAR = []; ieegChan = [];
    ieegCARAll = []; ieegChanAll = [];
    channels = [];
     specRaw = [];
     ieegSplit = [];
     ieegBase = [];
    [ieegSplit,~,trigOnset]=trialIEEGUpdate(Trials,1:length(allChannels),'ResponseStart','ieeg',tw.*1000);
    ieegSplit = permute(ieegSplit,[2,1,3]);
%     ieegBase = squeeze(trialIEEGUpdate(Trials,1:length(allChannels),'Start','ieeg',tw.*1000));
%     ieegBase = permute(ieegBase,[2,1,3]);
%     respId = find(~isnan(trigOnset));
%     ieegSplitResp = ieegSplit(:,respId,:);
%     ieegBaseResp = ieegBase(:,respId,:);
    
    
%% Bad channel removal
    ieegR=zeros(size(ieegSplit,1),size(ieegSplit,2)*size(ieegSplit,3));
for iChan=1:size(ieegSplit,1);
    ieegR(iChan,:)=reshape(ieegSplit(iChan,:,:),1,size(ieegSplit,2)*size(ieegSplit,3));
end
%ieegR=reshape(ieeg,size(ieeg,1),size(ieeg,2)*size(ieeg,3));
ieegR2=detrend(ieegR').^2;
iiZero=find(ieegR2==0);
ieegR2(iiZero)=.000000001;

ieegSTD=std(log(ieegR2),[],1);
ieegSTD=std(ieegR2,[],1);

[m s]=normfit(ieegSTD);
iiOutPlus1=find(ieegSTD>(3*s+m));
chanIn=setdiff(1:size(ieegSTD,2),iiOutPlus1);
[m s]=normfit(ieegSTD(chanIn));
iiOutPlus2=find(ieegSTD(chanIn)>(3*s+m));

badChannels=sort(cat(2,iiOutPlus1,chanIn(iiOutPlus2)));
%% Common average referencing
% ieegCarBase = carFilterImpedance(ieegBaseResp,badChannels);
% ieegCarSplit = carFilterImpedance(ieegSplitResp,badChannels);


[~,goodtrials] = remove_bad_trials(ieegSplit,8);
% goodTrialsCommon = extractCommonTrials(goodtrials);
% ieegBaseClean = ieegCarBase(:,goodTrialsCommon,:);
% ieegCarClean = ieegCarSplit(:,goodTrialsCommon,:);
 %% High Gamma Extraction 
 

fsDown =200;
specSamp = 80;
specDist = 50;
gInterval = 70;
normFactor = [];
ieegGammaPowerNorm= [];
[ieegGammaNorm,ieegGamma,p_masked] = ExtractHighGammaWrap(ieegSplit,fsD,fsDown,tw,etw,prtw,pstw,2);
timeGamma = linspace(etw(1),etw(2),size(ieegGammaNorm,3));
for iChan = 1:size(ieegGammaNorm,1)
       normFactor(iGamma,iChan,:) = [mean2(squeeze(ieegGammaBasedown(iChan,goodtrials{iChan},:))) std2(squeeze(ieegGammaBasedown(iGamma,goodtrials{iChan},:)))];
        %pChan(iChan) = permtest_sk((ieegGammaRespPower(iChan,goodtrials{iChan})),(ieegGammaBasePower(iChan,goodtrials{iChan})),10000);
        ieegGammaPowerNorm(iChan) = 20.*log10(mean2(squeeze(ieegGammaNorm(iChan,goodtrials{iChan},timeGamma>=pstw(1) & timeGamma<=pstw(2)))));
       % evokeSnr(iChan) = esnr(squeeze(ieegGammaRespdown(iChan,goodtrials{iChan},:)),squeeze(ieegGammaBasedown(iChan,goodtrials{iChan},:)));
end
   
%save(strcat(dLabels(iSubject).name,'_phonDecodeHGPower.mat'),'allChannels','p_masked','ieegGammaPowerNorm','motorChan', 'ifgChan');

% ieegGamma = [];
% 
% for iTrial = 1:size(ieegSplit,2)
%     iTrial
%     for iGamma = 1:length(gInterval)
%       [~,ieegTmp] = EcogExtractHighGammaTrial(double(squeeze(ieegSplitResp(:,iTrial,:))),fsD,fsDown,[gInterval(iGamma) gInterval(iGamma)+specSamp],tw,etwG,squeeze(normFactor(iGamma,:,:)),1);
%     ieegGamma(:,iTrial,iGamma,:) = ieegTmp;    
%     end
% end
% timeGamma = linspace(etwG(1),etwG(2),size(ieegGamma,4));


%% Channel Saving
sigChannel = [find(anatChan) ];
if(~isempty(sigChannel))
  [specRaw,pPerc,pProd]= getSpectrograms(ieegSplit,goodtrials,tw,etw,[1 200],prtw,pstw,[2.5 3],[70 150],fsD,1);
    [~, pvalsMCleanPerc] = fdr(pPerc,0.05);
    [ir,ic] = numSubplots(length(specRaw));
    figure;
    specPower = [];
for iChan = 1:length(specRaw)
     subplot(ir(1),ir(2),iChan);
% specChanMap(spec,[],1:length(channels),etw,[0.7 1.5],iChan)
% 
    [~,specPowerTemp] = specChanMap(specRaw,[],1:length(allChannels),find(pvalsMCleanPerc),etw,prtw,pstw,[1 200],[70 150],[0.7 1.4],iChan);
colormap(parula(4096));
title(allChannels(iChan));
specPower = [specPower; specPowerTemp];
 end  
  save(strcat(dLabels(iSubject).name,'_motorHGPack.mat'),'allChannels','p_masked','ieegGammaPowerNorm','anatChan','specPower','pvalsMCleanPerc');
end
end