addpath(genpath('E:\Box Sync\ECoG_Recon\matlab_code\'));
global DUKEDIR
DUKEDIR = 'E:\microECoG surgeries\S14 surgery';
dLabels = dir(DUKEDIR);
dLabels = dLabels(3:end);

tw = [-2 2]; % time window
etw = [-1.5 1.5];
etwG = [-1 1];% epoch time window
prtw = [-1.9 0.1]; % preonset time window
pstw = [0.25 0.75]; % postonset time window
gammaF = [70 150]; % frequency in Hz
load('channelMap.mat');
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
phon1 = phon1(respId);
phon2 = phon2(respId);
phon3 = phon3(respId);


 ieegSplitAll = cat(2,ieegSplit1);
 phonAll = [phon1 ];

  
trigOns = trigOnset(respId)./fsD;
trigOns2 = trigOnset2(respId)./fsD;
trigOfs = trigOffset(respId)./fsD;
%% Common average referencing
higImpId = find(log10(impedance1)>6);
ieegCarImp = carFilterImpedance(ieegSplitAll,higImpId);
ieegCarBase = carFilterImpedance(ieegBase,higImpId);
ieegCarAud = carFilterImpedance(ieegAuditory,higImpId);
%% High Gamma
load('sigChannel.mat');
fsDown = 200;
[~,goodtrials] = remove_bad_trials(double(ieegCarImp),14);
goodTrialsCommon = extractCommonTrials(goodtrials(pvalsMCleanProdResp));
ieegCarClean = ieegCarImp(:,goodTrialsCommon,:);
ieegCarAudClean = ieegCarAud(:,goodTrialsCommon,:);
ieegThetaBasedown = [];
ieegAlphaBasedown = [];
ieegBetaBasedown = [];
ieegGammaBasedown = [];
ieegHiGammaBasedown = [];
ieegHiGammaBaseNorm = [];
ieegGammaAllBasedown = [];
% for iTrial = 1:size(ieegCarAudClean,2)        
%         [~,ieegThetaBasedown(:,iTrial,:)] = EcogExtractHighGammaTrial(double(squeeze(ieegCarAudClean(:,iTrial,:))),fsD,fsDown,[4 7],tw,prtw,[]);
% end
% for iChan = 1:size(ieegThetaBasedown,1)
%     normFactor1(iChan,:) = [mean2(squeeze(ieegThetaBasedown(iChan,:,:))) std2(squeeze(ieegThetaBasedown(iChan,:,:)))];
%    % pChan(iChan) = permtest_sk(ieegGammaRespPower(iChan,:),ieegGammaBasePower(iChan,:),10000);
% end
% for iTrial = 1:size(ieegCarAudClean,2)        
%         [~,ieegAlphaBasedown(:,iTrial,:)] = EcogExtractHighGammaTrial(double(squeeze(ieegCarAudClean(:,iTrial,:))),fsD,fsDown,[7 14],tw,prtw,[]);
% end
% for iChan = 1:size(ieegAlphaBasedown,1)
%     normFactor2(iChan,:) = [mean2(squeeze(ieegAlphaBasedown(iChan,:,:))) std2(squeeze(ieegAlphaBasedown(iChan,:,:)))];
%    % pChan(iChan) = permtest_sk(ieegGammaRespPower(iChan,:),ieegGammaBasePower(iChan,:),10000);
% end
% for iTrial = 1:size(ieegCarAudClean,2)        
%         [~,ieegBetaBasedown(:,iTrial,:)] = EcogExtractHighGammaTrial(double(squeeze(ieegCarAudClean(:,iTrial,:))),fsD,fsDown,[15 30],tw,prtw,[]);
% end
% for iChan = 1:size(ieegBetaBasedown,1)
%     normFactor3(iChan,:) = [mean2(squeeze(ieegBetaBasedown(iChan,:,:))) std2(squeeze(ieegBetaBasedown(iChan,:,:)))];
%    % pChan(iChan) = permtest_sk(ieegGammaRespPower(iChan,:),ieegGammaBasePower(iChan,:),10000);
% end
% for iTrial = 1:size(ieegCarAudClean,2)        
%         [~,ieegGammaBasedown(:,iTrial,:)] = EcogExtractHighGammaTrial(double(squeeze(ieegCarAudClean(:,iTrial,:))),fsD,fsDown,[30 60],tw,prtw,[]);
% end
% for iChan = 1:size(ieegGammaBasedown,1)
%     normFactor4(iChan,:) = [mean2(squeeze(ieegGammaBasedown(iChan,:,:))) std2(squeeze(ieegGammaBasedown(iChan,:,:)))];
%    % pChan(iChan) = permtest_sk(ieegGammaRespPower(iChan,:),ieegGammaBasePower(iChan,:),10000);
% end
normFactorAll = [];
normFactor = [];
for iTrial = 1:size(ieegCarAudClean,2)
    
        [~,ieegHiGammaBasedown(:,iTrial,:),ieegGammaAllBasedown(:,iTrial,:,:)] = EcogExtractHighGammaTrial(double(squeeze(ieegCarAudClean(:,iTrial,:))),fsD,fsDown,[70 150],tw,prtw,[]);
end
for iChan = 1:size(ieegCarAudClean,1)
    iChan
    normFactor(iChan,:) = [mean2(squeeze(ieegHiGammaBasedown(iChan,:,:))) std2(squeeze(ieegHiGammaBasedown(iChan,:,:)))];
   % pChan(iChan) = permtest_sk(ieegGammaRespPower(iChan,:),ieegGammaBasePower(iChan,:),10000);
   for iFreq = 1:size(ieegGammaAllBasedown,4)
       normFactorAll(iChan,iFreq,:) = [mean2(squeeze(ieegGammaAllBasedown(iChan,:,:,iFreq))) std2(squeeze(ieegGammaAllBasedown(iChan,:,:,iFreq)))];
   end
end
% for iTrial = 1:size(ieegCarAudClean,2)        
%         [~,ieegHiGammaBaseNorm(:,iTrial,:)] = EcogExtractHighGammaTrial(double(squeeze(ieegCarAudClean(:,iTrial,:))),fsD,fsDown,[70 150],tw,prtw,normFactor5,1);
% end
%normFactor6 = [mean(ieegHiFreqBasedown(:,timeBaseDown>=5&timeBaseDown<=35),2) std(ieegHiFreqBasedown(:,timeBaseDown>=5&timeBaseDown<=35),0,2)];
ieegTheta = [];
ieegAlpha = [];
ieegBeta = [];
ieegGamma = [];
ieegHiGamma = [];
ieegGammaAll = [];
for iTrial = 1:size(ieegCarClean,2)
    iTrial
%     [~,ieegTmp] = EcogExtractHighGammaTrial(double(squeeze(ieegCarClean(:,iTrial,:))),fsD,fsDown,[4 7],tw,etwG,normFactor1,1);
%     ieegTheta(:,iTrial,:) = ieegTmp;
%     [~,ieegTmp] = EcogExtractHighGammaTrial(double(squeeze(ieegCarClean(:,iTrial,:))),fsD,fsDown,[7 14],tw,etwG,normFactor2,1);
%     ieegAlpha(:,iTrial,:) = ieegTmp;
%     [~,ieegTmp] = EcogExtractHighGammaTrial(double(squeeze(ieegCarClean(:,iTrial,:))),fsD,fsDown,[15 30],tw,etwG,normFactor3,1);
%     ieegBeta(:,iTrial,:) = ieegTmp;
%     [~,ieegTmp] = EcogExtractHighGammaTrial(double(squeeze(ieegCarClean(:,iTrial,:))),fsD,fsDown,[30 60],tw,etwG,normFactor4,1);
%     ieegGamma(:,iTrial,:) = ieegTmp;
    [~,ieegTmp,ieegTemp3] = EcogExtractHighGammaTrial(double(squeeze(ieegCarClean(:,iTrial,:))),fsD,fsDown,[70 150],tw,etwG,normFactor,1,normFactorAll);
    ieegHiGamma(:,iTrial,:) = ieegTmp;
    ieegGammaAll(:,iTrial,:,:) = ieegTemp3;
%     [~,ieegTmp] = EcogExtractHighGammaTrial(double(squeeze(ieegCarClean(:,iTrial,:))),fsD,fsDown,[185 250],tw,etwG,normFactor6,1);
%     ieegHiFreq(:,iTrial,:) = ieegTmp;
end
timeGamma = linspace(etwG(1),etwG(2),size(ieegHiGamma,3));
%% Phoneme Encoding
phonClean = phonAll(goodTrialsCommon);
binClass = [];
phonConsClass = [];
phonCVClass = [];
phonCvClassEnd = [];
phonIndClass = [];
phonIndClassEnd = [];
vowel = "aeiou";
for iTrial = 1:size(phonClean,2)
    [binClass,phonCVClass,phonIndClassTemp] = phonemeEncoder(phonClean{1,iTrial});
     phonIndClass(iTrial) = phonIndClassTemp;   
end
%%
timeSelect = timeGamma>=-0.5&timeGamma<=0.5;
ieegGammaSelect =(ieegHiGamma(:,:,timeSelect));
ieegGammaSelectAll =(ieegGammaAll(:,:,timeSelect,:));
% ieegGammaImage = [];
% for iGamma = 1:size(ieegGammaSelect,3)
%     ieegGammaImage(:,:,:,:,iGamma) = getImages(squeeze(ieegGammaSelect(:,:,iGamma,:)),chanMap,sigChannel);
% end

%save('gammaPhonemeRnn_1_5_seconds_v2.mat','ieegGammaSeries','ieegGammaImage','phonIndClass');
matSize = size(ieegGammaSelectAll);
 ieegModel = reshape(ieegGammaSelectAll,[matSize(1) matSize(2) matSize(3)*matSize(4)]);
            
labels = {'a','ae','i','u','b','p','v','g','k'};
accIter = [];
CmatCat = zeros(9,9)
for iTer = 1:5
[~,ytestAll,ypredAll] = stmfDecodeWrap(ieegGammaSelect,phonIndClass,[0 1],[0 1],20,0);
%[~,ytestAll,ypredAll,~,stmfTemplate] = stmfDecodeWrap(ieegModel,phonIndClass,[0 1],[0 1],20,0);
CmatAll = confusionmat(ytestAll,ypredAll);
acc = trace(CmatAll)/sum(sum(CmatAll))
accIter(iTer) = acc;
CmatCat = CmatCat+CmatAll;
% stmfmean = squeeze(mean(stmfTemplate,1));
% stmfmean = reshape(stmfmean,[size(stmfmean,1) matSize(1) matSize(3)*matSize(4)]);
% stmfSpace = squeeze(mean(stmfmean,3));
end
% for iPhone = 1:9
% figure;
% chanView(stmfSpace(iPhone,:),chanMap,selectedChannels,isnan(chanMap),labels(iPhone),[],[],0);
% end
accAll = trace(CmatCat)/sum(sum(CmatCat))
figure; 
cm = confusionchart(CmatCat,labels,'Normalization','row-normalized');
sortClasses(cm,labels);
set(gca,'FontSize',20);
title(['Accuracy : ' num2str(accAll*100) '%']);

[phonError,cmatvect,phonemeDistVect] = phonemeDistanceError(CmatCat);
cmatvect(phonemeDistVect==0) = [];
phonemeDistVect(phonemeDistVect==0) = [];

distMod = fitlm(phonemeDistVect,cmatvect)
figure;
scatter(phonemeDistVect,cmatvect,'filled')
hold on;
plot(phonemeDistVect,predict(distMod,phonemeDistVect));
xlabel('Phoneme distance (bits)')
ylabel('Decoding error');
%xlim([0.5 9.5]);
set(gca,'FontSize',15);

%%

