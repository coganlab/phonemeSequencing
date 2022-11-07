addpath(genpath('E:\Box Sync\ECoG_Recon\matlab_code\'));
global DUKEDIR
DUKEDIR = 'E:\microECoG surgeries\S14 surgery';
dLabels = dir(DUKEDIR);
dLabels = dLabels(3:end);

tw = [-2 2.5]; % time window
etw = [-1.5 1.5];
etwG = [-1 1];% epoch time window
etwGAuditory = [0 2];
prtw = [-0.5 0]; % preonset time window
pstw = [0.25 0.75]; % postonset time window
gammaF = [70 150]; % frequency in Hz
load('channelMap.mat');
chanMap = fliplr(chanMap);
selectedChannels = sort(chanMap(~isnan(chanMap)))';
%% data Loading
Experiment = loadExperiment('S14');
fsD = Experiment.processing.ieeg.sample_rate;
Trials = dbTrials('S14',Experiment.recording.recording_day,'Speech_OvertMimeMove');

trialFiles = strcat('\',Experiment.recording.recording_day,'\mat\trialInfo.mat');
[ieegSplit1,~,trigOnset]=trialIEEGUpdate(Trials,selectedChannels,'phon1Onset','ieeg',tw.*1000);
[ieegSplit2,~,trigOnset2]=trialIEEGUpdate(Trials,selectedChannels,'phon2Onset','ieeg',tw.*1000);
[ieegSplit3,~,~]=trialIEEGUpdate(Trials,selectedChannels,'phon3Onset','ieeg',tw.*1000);
[~,~,trigOffset]=trialIEEGUpdate(Trials,selectedChannels,'ResponseOffset','ieeg',tw.*1000);
[ieegAuditory,~,trigAuditory]=trialIEEGUpdate(Trials,selectedChannels,'Auditory','ieeg',tw.*1000);
ieegAuditory = permute(ieegAuditory,[2,1,3]);

ieegSplit1 = permute(ieegSplit1,[2,1,3]);
ieegSplit2 = permute(ieegSplit2,[2,1,3]);
ieegSplit3 = permute(ieegSplit3,[2,1,3]);
phon1 = [Trials.phon1];
phon2 = [Trials.phon2];
phon3 = [Trials.phon3];
ieegBase = squeeze(trialIEEGUpdate(Trials(end),selectedChannels,'ResponseOnset','ieeg',[2 42].*1000));
 
respId = find(~isnan(trigOnset));
ieegSplit1 = ieegSplit1(:,respId,:);
ieegSplit2 = ieegSplit2(:,respId,:);
ieegSplit3 = ieegSplit3(:,respId,:);
ieegAuditory = ieegAuditory(:,respId,:);
phon1 = phon1(respId);
phon2 = phon2(respId);
phon3 = phon3(respId);


 ieegSplitAll = cat(2,ieegSplit1);
 phonAll = [phon1; phon2; phon3 ];
phonSequence =strcat(phon1,phon2,phon3);
  
trigOns = trigOnset(respId)./fsD;
trigOns2 = trigOnset2(respId)./fsD;
trigOfs = trigOffset(respId)./fsD;
%% Common average referencing
load('impedance.mat')
higImpId = find(log10(impedance1)>6);
ieegCarImp = carFilterImpedance(ieegSplitAll,higImpId);
ieegCarBase = carFilterImpedance(ieegBase,higImpId);
ieegCarAud = carFilterImpedance(ieegAuditory,higImpId);
% ieegCarBase = filt60(double(ieegCarBase),fsD);
% ieegCarImpFilt = [];
% for iChan = 1:size(ieegCarImp,1)
%     ieegCarImpFilt(iChan,:,:) = filt60(squeeze(double(ieegCarImp(iChan,:,:))),fsD);
% end
% ieegCarImp = ieegCarImpFilt;
%% Spectrogram 
[~,goodtrialsbase] = remove_bad_trials(ieegCarAud,14);
[specCarBase]= getSpectrograms(ieegCarAud,goodtrialsbase,tw,etw,[1 200],prtw,pstw,[70 150],fsD,1);
 
figure;
[~,~,meanFreqChanOut] = specChanMap(specCarBase,(chanMap),selectedChannels,[],etw,prtw,pstw,[1 200],[70 150],[0.7 1.4],0,[]);
[~,goodtrials] = remove_bad_trials(ieegCarImp,14);
[specCarResp,pvalCarProd]= getSpectrograms(ieegCarImp,goodtrials,tw,etw,[1 200],[-1.5 -1],[-0.25 0.25],[70 150],fsD,1);
[~, pvalsCarProd] = fdr(pvalCarProd,0.05);
%     [~, pvalsMCleanProd] = fdr(pProd,0.05);
figure;
specChanMap(specCarResp,fliplr(chanMap),selectedChannels,find(pvalsCarProd),etw,[-1.5 -1],[-0.25 0.25],[1 200],[70 150],[0.7 1.4],0,meanFreqChanOut);





%% High Gamma Filtering
load('S14_sigChannel.mat');

[~,goodtrials] = remove_bad_trials(ieegCarImp,14);
goodTrialsCommon = extractCommonTrials(goodtrials(sigChannel));
ieegCarClean = ieegCarImp(:,goodTrialsCommon,:);
ieegCarAudClean = ieegCarAud(:,goodTrialsCommon,:);
fsDown =200;

gInterval = [70 150];
normFactor = [];

    ieegGammaBasedown = [];
    ieegGammaRespdown = [];
    for iTrial = 1:size(ieegCarAudClean,2)
        iTrial
        [~,ieegGammaBasedown(:,iTrial,:)] = EcogExtractHighGammaTrial(double(squeeze(ieegCarAudClean(:,iTrial,:))),fsD,fsDown,[gInterval(1) gInterval(2)],tw,prtw,[]);
        %[~,ieegGammaRespdown(:,iTrial,:)] = EcogExtractHighGammaTrial(double(squeeze(ieegCarClean(:,iTrial,:))),fsD,fsDown,[gInterval(iGamma) gInterval(iGamma)+specSamp],tw,[-0.25 0.25],[]);
    end
%     ieegGammaRespPower = squeeze(mean(ieegGammaRespdown,3));
%     ieegGammaBasePower = squeeze(mean(ieegGammaBasedown,3));
    for iChan = 1:size(ieegGammaBasedown,1)
    normFactor(iChan,:) = [mean2(squeeze(ieegGammaBasedown(iChan,:,:))) std2(squeeze(ieegGammaBasedown(iChan,:,:)))];
   % pChan(iChan) = permtest_sk(ieegGammaRespPower(iChan,:),ieegGammaBasePower(iChan,:),10000);
    end
     %[p_fdr, p_masked] = fdr( pChan, 0.05);


ieegGammaClean = [];
ieegGammaAud = [];
for iTrial = 1:size(ieegCarClean,2)
    iTrial
    [~,ieegTmp] = EcogExtractHighGammaTrial(double(squeeze(ieegCarAudClean(:,iTrial,:))),fsD,fsDown,[gInterval(1) gInterval(2)],tw,etwGAuditory,squeeze(normFactor),2);
        ieegGammaAud(:,iTrial,:) = ieegTmp;
    
        [~,ieegTmp] = EcogExtractHighGammaTrial(double(squeeze(ieegCarClean(:,iTrial,:))),fsD,fsDown,[gInterval(1) gInterval(2)],tw,etwG,squeeze(normFactor),2);
        ieegGammaClean(:,iTrial,:) = ieegTmp;    
    
end
timeGamma = linspace(etwG(1),etwG(2),size(ieegGammaClean,3));



%travellingWaveMovie(ieegGammaRespMean,fsD,chanMap,selectedChannels,timeGamma,[-1 0.5],[0 2],240,'S14GammaResponse','z-score');

%% Phoneme Encoding
phonClean = phonAll(:,goodTrialsCommon);
phonSequenceClean = phonSequence(:,goodTrialsCommon);
binClass = [];
phonConsClass = [];
phonCVClass = [];
phonCvClassEnd = [];
phonIndClass = [];
phonIndClassEnd = [];
phontactic =[];
vowel = "aeiou";
for iTrial = 1:size(phonClean,2)
    [binClassTemp,phonCVClassTemp,phonIndClassTemp] = phonemeEncoder(phonClean{1,iTrial});
    phonIndClass(iTrial,1) = phonIndClassTemp;     
    phonCVClass(iTrial,1) = phonCVClassTemp;
    [binClassTemp,phonCVClassTemp,phonIndClassTemp] = phonemeEncoder(phonClean{2,iTrial});
    phonIndClass(iTrial,2) = phonIndClassTemp;
    phonCVClass(iTrial,2) = phonCVClassTemp;
    [binClassTemp,phonCVClassTemp,phonIndClassTemp] = phonemeEncoder(phonClean{3,iTrial});
    phonIndClass(iTrial,3) = phonIndClassTemp;
    phonCVClass(iTrial,3) = phonCVClassTemp;
    phonid = find(strcmp([PhonemeSequencingInfoS1.TOKEN],phonSequenceClean{iTrial}));
     
    if(isempty(phonid))
        phonSequenceClean{iTrial}
        phontactic(iTrial,:) = nan(1,size(PhonemeSequencingInfoS1,10));
    else
        phontactic(iTrial,:) = table2array(PhonemeSequencingInfoS1(phonid,2:11));
    end
    
    
end

        %% Individual phoneme decoding
  phontactic(phontactic==0) = nan;
        phontacticLog = -log(phontactic(:,1:9));
        phontacticLog(:,10) = phontactic(:,10);
cvcIds = phonIndClass(:,1)>4;
 for tRange = 0.5
timeSelect = timeGamma>=-tRange&timeGamma<=tRange;
ieegGammaSelect =(ieegGammaClean(featId==2,:,timeSelect));

accIter = []; phonError = [];
phonTacVariables = PhonemeSequencingInfoS1.Properties.VariableNames(2:11);
for iPhon = 1:10
    phonVariable = phontacticLog(:,iPhon)';
    for iTer = 1
        iTer
        [ytestAll,ypredAll,optimDimAll] = pcaLinearRegressDecoderWrap(ieegGammaSelect,phonVariable,[0 1],[0 1],80,20);
         %[~,ytestAll,ypredAll,~,stmfTemplate] = stmfDecodeWrap(ieegModel,phonIndClass(:,1)',[0 1],[0 1],20,0);
        %[~,ytestAll,ypredAll] = linearDecoder(ieegModel,phonIndClass(:,1)',[0 1],[0 1],10,0);
        rsqrd = 1 - nansum((ytestAll - ypredAll).^2)/nansum((ytestAll - nanmean(ytestAll)).^2);
        distMod = fitlm(ytestAll,ypredAll)
        figure;
        scatter(ytestAll,ypredAll,'filled')
        hold on;
        plot(ytestAll,predict(distMod,ytestAll'));
        xlabel('True Phonotactics ')
        ylabel('Predicted Phonotactics');
        %xlim([0.5 9.5]);
        set(gca,'FontSize',15);
        title([phonTacVariables(iPhon) ' R-squared : ' num2str(rsqrd)]);
    end
end
 end

 
%% Time-series prediction
timeRange = -1:0.01:0.75;
phontactic(phontactic==0) = nan;
phontacticLog = -log(phontactic(:,1:9));
phontacticLog(:,10) = phontactic(:,10);
accTime = zeros(1,length(timeRange));
pValTime = zeros(1,length(timeRange));
aucTime = [];
%ieegGammaSelect = squeeze(ieegGamma(sigChannel,goodTrialsCommon,1,:));
for iPhon = 1:10
    for iTimeTrain = 1:length(timeRange)
        iTimeTrain
        tsTrain = [timeRange(iTimeTrain) timeRange(iTimeTrain)+0.25];
            [ytestAll,ypredAll,optimDimAll] = pcaLinearRegressDecoderWrap(ieegGammaClean(sigChannel,cvcIds,:),(phontacticLog(cvcIds,iPhon)'),etwG,tsTrain,[80],20);
        distMod = fitlm(ytestAll,ypredAll);
        accTime(iPhon,iTimeTrain) = distMod.Rsquared.Ordinary;
        pValTime(iPhon,iTimeTrain) = distMod.Coefficients.pValue(2);
    end
end


timeRangeCorrect = timeRange+0.125;
%          figure;
%          plot(timeRangeCorrect,accTime);
%          hold on;
%          [pvalnew,p_masked] = fdr(pValTime(:),0.01);
%          scatter(timeRangeCorrect(pValTime<pvalnew),0.2.*ones(1,sum(pValTime<pvalnew)),'filled');
%          axis square;
%          set(gca,'FontSize',20);
%          xlabel('Time from production onset (s)');
%          ylabel('BLICK prediction R^2');
 
 
 figure; 
 for iPhon = 1:10
     subplot(2,5,iPhon);
     plot(timeRangeCorrect,accTime(iPhon,:),'LineWidth',2);
     hold on;
     [pvalnew,p_masked] = fdr(pValTime(iPhon,:),0.01);
         scatter(timeRangeCorrect(pValTime(iPhon,:)<pvalnew),0.5.*ones(1,sum(pValTime(iPhon,:)<pvalnew)),'filled');
     xlabel('Time (s)')
     title(phonTacVariables(iPhon));
     set(gca,'FontSize',15);
 end
 
 
  %% Individual phoneme decoding - time series
  
  timeRange = -1.5:0.01:1.5;

accTime = zeros(10,length(timeRange));
  
 for iTime = 1:length(timeRange)
timeSelect = timeGamma>=timeRange(iTime)&timeGamma<=timeRange(iTime)+0.25;
ieegGammaSelect =(ieegGammaClean(:,:,1,timeSelect));
% ieegGammaImage = [];
% for iGamma = 1:size(ieegGammaSelect,3)
%     ieegGammaImage(:,:,:,:,iGamma) = getImages(squeeze(ieegGammaSelect(:,:,iGamma,:)),chanMap,find(pvalsMCleanProdResp));
% end
% ieegGammaSeries = permute(ieegGammaSelect,[2 4 1 3]);
% save('gammaPhonemeRnn_1_seconds_v2.mat','ieegGammaSeries','ieegGammaImage','phonIndClass');
matSize = size(ieegGammaSelect);
ieegModel = reshape(ieegGammaSelect,[matSize(1) matSize(2) matSize(3)*matSize(4)]);
%             

accIter = []; phonError = [];
phonTacVariables = PhonemeSequencingInfoS1.Properties.VariableNames(2:11);
for iPhon = 1:10
    iPhon
    phonVariable = phontactic(:,iPhon)';
    for iTer = 1
        
        [ytestAll,ypredAll,optimDimAll] = pcaLinearRegressDecoderWrap(ieegModel,phonVariable,[0 1],[0 1],80,20);
         %[~,ytestAll,ypredAll,~,stmfTemplate] = stmfDecodeWrap(ieegModel,phonIndClass(:,1)',[0 1],[0 1],20,0);
        %[~,ytestAll,ypredAll] = linearDecoder(ieegModel,phonIndClass(:,1)',[0 1],[0 1],10,0);
        distMod = fitlm(ytestAll,ypredAll);
        accTime(iPhon,iTime) = 1 - nansum((ytestAll - ypredAll).^2)/nansum((ytestAll - nanmean(ytestAll)).^2);
%         figure;
%         scatter(ytestAll,ypredAll,'filled')
%         hold on;
%         plot(ytestAll,predict(distMod,ytestAll'));
%         xlabel('True Phonotactics ')
%         ylabel('Predicted Phonotactics');
%         %xlim([0.5 9.5]);
%         set(gca,'FontSize',15);
%         title([phonTacVariables(iPhon) ' R-squared : ' num2str(distMod.Rsquared.Ordinary)]);
    end
end
 end
 
 
 figure; 
 for iPhon = 1:9
     subplot(3,3,iPhon);
     plot(timeRange,accTime(iPhon,:),'LineWidth',2);
     xlabel('Time (s)')
     title(phonTacVariables(iPhon));
     set(gca,'FontSize',15);
 end
 