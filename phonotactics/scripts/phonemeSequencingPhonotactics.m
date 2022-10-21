addpath(genpath('E:\Box Sync\Box Sync\ECoG_Recon\matlab_code\'));
global DUKEDIR
DUKEDIR = 'E:\Box Sync\Box Sync\CoganLab\D_Data\Phoneme_Sequencing';
dLabels = dir(DUKEDIR);
dLabels = dLabels(3:end);

tw = [-2 2]; % time window
etwG = [-1 1];
etwGAuditory = [0 2];% epoch time window
prtw = [-0.5 0]; % preonset time window
pstw = [0.25 0.75]; % postonset time window
gammaF = [70 150]; % frequency in Hz
%%
for iSubject =9
    
   d = []; ieegCarResp = []; ieegCarImpFilt = []; ieegGammaResp = []; ieegResponse = [];
    Experiment = loadExperiment(dLabels(iSubject).name);
    fsD = Experiment.recording.sample_rate;
    Trials = dbTrials(dLabels(iSubject).name,Experiment.recording.recording_day,'Speech_OvertMimeMove');
    allChannels = string({Experiment.channels.name});
    for iChan = 1:length(allChannels)
        allChannels(iChan) = strcat(dLabels(iSubject).name,'-',allChannels(iChan));
    end
    trialFiles = strcat('\',Experiment.recording.recording_day,'\mat\trialInfo.mat');
    load([DUKEDIR '\' dLabels(iSubject).name '\' trialFiles])
    [~,chan2select] = intersect(allChannels,prodDelayElecs);
    if(isempty(chan2select))
        continue;
    end
    

%% Loading all the data
    
 disp('Loading data...')   
     
    [ieegSplit,~,trigOnset]=trialIEEGUpdate(Trials,1:length(allChannels),'ResponseStart','ieeg',tw.*1000);
    ieegSplit = permute(ieegSplit,[2,1,3]);
    ieegBase = squeeze(trialIEEGUpdate(Trials,1:length(allChannels),'Start','ieeg',tw.*1000));
    ieegBase = permute(ieegBase,[2,1,3]);
    ieegAuditory = squeeze(trialIEEGUpdate(Trials,1:length(allChannels),'Auditory','ieeg',tw.*1000));
    ieegAuditory = permute(ieegAuditory,[2,1,3]);
    respId = find(~isnan(trigOnset));
    ieegSplitResp = ieegSplit(:,respId,:);
    ieegBaseResp = ieegBase(:,respId,:);
    ieegAudResp = ieegAuditory(:,respId,:);
  disp('done');  
%% Bad channel removal
disp('removing bad channels')
    ieegR=zeros(size(ieegSplitResp,1),size(ieegSplitResp,2)*size(ieegSplitResp,3));
for iChan=1:size(ieegSplitResp,1);
    ieegR(iChan,:)=reshape(ieegSplitResp(iChan,:,:),1,size(ieegSplitResp,2)*size(ieegSplitResp,3));
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
disp('done')
%% Common average referencing
disp('common average referencing')
ieegCarBase = carFilterImpedance(ieegBaseResp,badChannels);
ieegCarSplit = carFilterImpedance(ieegSplitResp,badChannels);
ieegCarAud = carFilterImpedance(ieegAudResp,badChannels);
% ieegCarBase = ieegBaseResp;
% ieegCarSplit = ieegSplitResp;
goodChannels = setdiff(1:size(ieegSplitResp,1),badChannels);
% ieegCarBase = ieegBaseResp;
% ieegCarSplit = ieegSplitResp;
[~,goodtrials] = remove_bad_trials(ieegCarSplit,14);
goodTrialsCommon = extractCommonTrials(goodtrials(chan2select));
ieegBaseClean = ieegCarBase(chan2select,goodTrialsCommon,:);
ieegCarClean = ieegCarSplit(chan2select,goodTrialsCommon,:);
ieegAudClean = ieegCarAud(chan2select,goodTrialsCommon,:);
disp('done')
 %% High Gamma Extraction 
 disp('Extracting High Gamma')
fsDown =200;
gInterval = [70 150];
normFactor = [];


    ieegGammaBasedown = [];
    ieegGammaRespdown = [];
    for iTrial = 1:size(ieegBaseClean,2)
        
        [~,ieegGammaBasedown(:,iTrial,:)] = EcogExtractHighGammaTrial(double(squeeze(ieegBaseClean(:,iTrial,:))),fsD,fsDown,[gInterval(1) gInterval(2)],tw,prtw,[],[]);
        % [~,ieegGammaRespdown(:,iTrial,:)] = EcogExtractHighGammaTrial(double(squeeze(ieegCarClean(:,iTrial,:))),fsD,fsDown,[gInterval(1) gInterval(2)],tw,[-0.25 0.25],[]);
    end
%     
%     ieegGammaRespPower = log10(squeeze(mean(ieegGammaRespdown.^2,3)));
%     ieegGammaBasePower = log10(squeeze(mean(ieegGammaBasedown.^2,3)));
    pChan = [];
%     ieegGammaPowerNorm = [];
%     evokeSnr = []
    for iChan = 1:size(ieegGammaBasedown,1)
        normFactor(iChan,:) = [mean2(squeeze(ieegGammaBasedown(iChan,:,:))) std2(squeeze(ieegGammaBasedown(iChan,:,:)))];
%         pChan(iChan) = permtest_sk((ieegGammaRespPower(iChan,:)),(ieegGammaBasePower(iChan,:)),10000);
%         ieegGammaPowerNorm(iChan) = 10.*log10(mean(ieegGammaRespPower(iChan,goodtrials{iChan})./mean(ieegGammaBasePower(iChan,goodtrials{iChan}))));
%         evokeSnr(iChan) = esnr(squeeze(ieegGammaRespdown(iChan,goodtrials{iChan},:)),squeeze(ieegGammaBasedown(iChan,goodtrials{iChan},:)));
    end
 %    [p_fdr, pvalsMCleanProd] = fdr( pChan, 0.01);

   
%save(strcat(dLabels(iSubject).name,'_phonDecodeHGPower.mat'),'allChannels','p_masked','ieegGammaPowerNorm','motorChan', 'ifgChan');

ieegGamma = [];
ieegGammaAud = [];
for iTrial = 1:size(ieegCarClean,2)
    
    [~,ieegTmp] = EcogExtractHighGammaTrial(double(squeeze(ieegAudClean(:,iTrial,:))),fsD,fsDown,[gInterval(1) gInterval(2)],tw,etwGAuditory,squeeze(normFactor),1);
        ieegGammaAud(:,iTrial,:) = ieegTmp;
      [~,ieegTmp] = EcogExtractHighGammaTrial(double(squeeze(ieegCarClean(:,iTrial,:))),fsD,fsDown,[gInterval(1) gInterval(2)],tw,etwG,squeeze(normFactor),1);
        ieegGamma(:,iTrial,:) = ieegTmp;    
    
end
ieegGamma = squeeze(ieegGamma);
timeGamma = linspace(etwG(1),etwG(2),size(ieegGamma,3));
disp('done');

%% Phoneme encoding
disp('phoneme coding')
trialNames = [];

trialInfoResp = trialInfo(respId);
trialClean = trialInfoResp(goodTrialsCommon);


phonPairAll = [];

binClass = [];
phonConsClass = [];
phonCVClass = [];
phonCvClassEnd = [];
phonIndClass = [];
phonIndClassEnd = [];
phontactic = [];
for iTrial = 1:length(trialClean)
    phonPairParse = [];
    if(iscell(trialClean))
        trialNames = (trialClean{iTrial}.sound(1:end-4));
    else
        trialNames = (trialClean(iTrial).sound(1:end-4));
    end
    trialNamesParse = strrep(trialNames,'ae','z');
    trialNamesParse = num2cell(trialNamesParse);
    trialNamesParse = strrep(trialNamesParse,'z','ae');
    for iPhon = 1:3
         [~,phonCVClass(iTrial,iPhon),phonIndClass(iTrial,iPhon)] = phonemeEncoder(trialNamesParse{iPhon});          
    end
     phonid = find(strcmp([PhonemeSequencingInfoS1.TOKEN],trialNames));
     
    if(isempty(phonid))
        trialNames
        phontactic(iTrial,:) = nan(1,size(PhonemeSequencingInfoS1,10));
    else
        phontactic(iTrial,:) = table2array(PhonemeSequencingInfoS1(phonid,2:11));
    end
end
disp('done')

%% CVC
disp('phonotactic modeling')
cvcIds = (phonIndClass(:,1)>4);
r2all = []; 
accVcv = [];
r2Chan = [];
accVcvChan = [];
phonError = [];
phonTacVariables = PhonemeSequencingInfoS1.Properties.VariableNames(2:11);
for iPhon = 10
    phonVariable = (phontactic(:,iPhon)');
    
        disp('Multivariate')
        [ytestAll,ypredAll,optimDimAll] = pcaLinearRegressDecoderWrap(ieegGammaAud,phonVariable,etwGAuditory,[0.5 1],[10:10:90],10);
         %[~,ytestAll,ypredAll,~,stmfTemplate] = stmfDecodeWrap(ieegModel,phonIndClass(:,1)',[0 1],[0 1],20,0);
        %[~,ytestAll,ypredAll] = linearDecoder(ieegModel,phonIndClass(:,1)',[0 1],[0 1],10,0);
        distMod = fitlm(ytestAll,ypredAll);
         r2all = distMod.Rsquared.Ordinary;
         disp('done')
%          disp('Univariate')
%          for iChan = 1:size(ieegGamma,1)
%              [ytestAll,ypredAll] = pcaLinearRegressDecoderWrap(ieegGamma(iChan,:,:),phonVariable,etwG,[-0.75 0.75],[10:10:90],10);
%              %[~,ytestAll,ypredAll,~,stmfTemplate] = stmfDecodeWrap(ieegModel,phonIndClass(:,1)',[0 1],[0 1],20,0);
%             %[~,ytestAll,ypredAll] = linearDecoder(ieegModel,phonIndClass(:,1)',[0 1],[0 1],10,0);
%             distMod = fitlm(ytestAll,ypredAll);
%              r2Chan(iChan) = distMod.Rsquared.Ordinary;
%          end
%          disp('done')
%         figure;
%         scatter(ytestAll,ypredAll,'filled')
%         hold on;
%         plot(ytestAll,predict(distMod,ytestAll'));
%         xlabel('True Phonotactics ')
%         ylabel('Predicted Phonotactics');
%         %xlim([0.5 9.5]);
%         set(gca,'FontSize',15);
%        
        
    
end
%% Time-series prediction
timeRange = 0:0.01:1.75;

accTime = zeros(1,length(timeRange));
pValTime = zeros(1,length(timeRange));
aucTime = [];
%ieegGammaSelect = squeeze(ieegGamma(sigChannel,goodTrialsCommon,1,:));
for iTimeTrain = 1:length(timeRange)
    iTimeTrain
    tsTrain = [timeRange(iTimeTrain) timeRange(iTimeTrain)+0.25];
        [ytestAll,ypredAll,optimDimAll] = pcaLinearRegressDecoderWrap(ieegGammaAud,phonVariable,etwGAuditory,tsTrain,[10:10:90],10);
    distMod = fitlm(ytestAll,ypredAll);
    accTime(iTimeTrain) = distMod.Rsquared.Ordinary;
    pValTime(iTimeTrain) = distMod.Coefficients.pValue(2);
end


timeRangeCorrect = timeRange;
         figure;
         plot(timeRangeCorrect,accTime);
         hold on;
         [pvalnew,p_masked] = fdr(pValTime(:),0.01);
         scatter(timeRangeCorrect(pValTime<pvalnew),0.2.*ones(1,sum(pValTime<pvalnew)),'filled');
         axis square;
         set(gca,'FontSize',20);
         xlabel('Time from Auditory onset (s)');
         ylabel('BLICK prediction R^2');
         
%%
disp('Auditory temporal generalization')
timeRange = 0:0.025:1.75;

accTime = zeros(length(timeRange),length(timeRange));
pValTime = zeros(length(timeRange),length(timeRange));
aucTime = [];
%ieegGammaSelect = squeeze(ieegGamma(sigChannel,goodTrialsCommon,1,:));
for iTimeTrain = 1:length(timeRange)
    
    tsTrain = [timeRange(iTimeTrain) timeRange(iTimeTrain)+0.25];
    for iTimeTest = 1:length(timeRange)        
        tsTest = [timeRange(iTimeTest) timeRange(iTimeTest)+0.25];    
        [ytestAll,ypredAll] = pcaLinearRegressDecoderWrapTrainTest(ieegGammaAud,phonVariable,etwGAuditory,tsTrain,tsTest,mean(optimDimAll),10);
        distMod = fitlm(ytestAll,ypredAll);
        accTime(iTimeTrain,iTimeTest) = distMod.Rsquared.Ordinary;
     
    end
         save(strcat(dLabels(iSubject).name,'_auditory_vs_auditory_v1.mat'),'accTime','timeRange');
end
disp('done')
%%
disp('Response temporal generalization')
timeRange = -1:0.025:0.75;

accTime = zeros(length(timeRange),length(timeRange));
pValTime = zeros(length(timeRange),length(timeRange));
aucTime = [];
%ieegGammaSelect = squeeze(ieegGamma(sigChannel,goodTrialsCommon,1,:));
for iTimeTrain = 1:length(timeRange)
    
    tsTrain = [timeRange(iTimeTrain) timeRange(iTimeTrain)+0.25];
    for iTimeTest = 1:length(timeRange)        
        tsTest = [timeRange(iTimeTest) timeRange(iTimeTest)+0.25];    
        [ytestAll,ypredAll] = pcaLinearRegressDecoderWrapTrainTest(ieegGamma,phonVariable,etwG,tsTrain,tsTest,mean(optimDimAll),10);
        distMod = fitlm(ytestAll,ypredAll);
        accTime(iTimeTrain,iTimeTest) = distMod.Rsquared.Ordinary;
     
    end
         save(strcat(dLabels(iSubject).name,'_Response_vs_Response_v1.mat'),'accTime','timeRange');
end
disp('done')
    %%
disp(strcat(dLabels(iSubject).name, ' saving'));
save(strcat(dLabels(iSubject).name,'_phonotactics_v2.mat'),'r2all','r2Chan');
end