addpath(genpath('E:\Box Sync\ECoG_Recon\matlab_code\'));
global DUKEDIR
DUKEDIR = 'C:\Users\sd355\Downloads\phoneme-sequence';
dLabels = dir(DUKEDIR);
dLabels = dLabels(3:end);

tw = [-2 2]; % time window
etw = [-1.5 1.5];
etwG = [-1 1];% epoch time window
prtw = [-0.5 0]; % preonset time window
pstw = [0.25 0.75]; % postonset time window
gammaF = [70 150]; % frequency in Hz
%%
for iSubject = 1
    d = []; ieegCarImp = []; ieegCarImpFilt = []; ieegGamma = []; ieegSplit = [];
    Experiment = loadExperiment(dLabels(iSubject).name);
    fsD = Experiment.recording.sample_rate;
    Trials = dbTrials(dLabels(iSubject).name,Experiment.recording.recording_day,'Speech_OvertMimeMove');
    allChannels = string({Experiment.channels.name});
    trialFiles = strcat('\',Experiment.recording.recording_day,'\mat\trialInfo.mat');
    load([DUKEDIR '\' dLabels(iSubject).name '\' trialFiles])
    %% Recon Visualization
%     iSubject = 1;
%     dLabels(iSubject).name = 'D2';
    cfg = [];
    % cfg.alpha = 0.4;
    % cfg.use_brainshifted = 1;
    handles = plot_subj(dLabels(iSubject).name, cfg);
    clabels= list_electrodes({dLabels(iSubject).name});
    inUnitFile = clabels';
    gridIndex = zeros(size(inUnitFile));
    %% Preprocessing
    ieegCAR = []; ieegChan = [];
    ieegCARAll = []; ieegChanAll = [];
    channels = [];
     specRaw = [];
     ieegSplit = [];
     ieegBase = [];
    [ieegSplit,~,trigOnset]=trialIEEGUpdate(Trials,1:length(allChannels),'ResponseStart',tw.*1000);
    ieegSplit = permute(ieegSplit,[2,1,3]);
    ieegBase = squeeze(trialIEEGUpdate(Trials(1),1:length(allChannels),'ResponseStart',[-40 0].*1000));
    
    respId = find(~isnan(trigOnset));
    ieegSplit = ieegSplit(:,respId,:);
    trigOns = trigOnset(respId);
    
    
    rmsVal = std(ieegBase,0,2);
    higRmsId = rmsVal>mean(rmsVal)+ 3*std(rmsVal);
    ieegCarBase = carFilterImpedance(ieegBase,higRmsId);
    ieegCarBase = filt60(double(ieegCarBase),fsD);
    %ieegCarImp = carFilterImpedance(ieegSplit,higRmsId);
    ieegCarImpFilt = ieegSplit;
    for iChan = 1:size(ieegSplit,1)
        ieegCarImpFilt(iChan,:,:) = filt60(squeeze(double(ieegSplit(iChan,:,:))),fsD);
    end
    ieegSplit = ieegCarImpFilt;
    load(strcat(dLabels(iSubject).name,'_sigChannel.mat'))
    time = linspace(tw(1),tw(2),size(ieegSplit,3));
    [~,goodtrials] = remove_bad_trials(ieegSplit,10);
    
   %load(strcat(dLabels(iSubject).name,'_sigChannel.mat'))
    [specRaw,pPerc,pProd]= getSpectrograms(ieegSplit,[],tw,[-1.5 1.5],[1 200],[1 1.5],[-0.25 0.25],[2.5 3],[70 150],fsD,1);
    [~, pvalsMCleanPerc] = fdr(pPerc,0.01);
    [ir,ic] = numSubplots(length(specRaw));
   % goodTrialsCommon = extractCommonTrials(goodtrials(pvalsMCleanPerc));
    figure;
    specPower = [];
for iChan = 1:length(specRaw)
     subplot(ir(1),ir(2),iChan);
% specChanMap(spec,[],1:length(channels),etw,[0.7 1.5],iChan)
% 
    [~,specPowerTemp] = specChanMap(specRaw,[],1:length(allChannels),find(pvalsMCleanPerc),[-1.5 1.5],[1 1.5],[-0.25 0.25],[1 200],[70 150],[0.7 1.4],iChan);
colormap(parula(4096));
title(allChannels(iChan));
specPower = [specPower; specPowerTemp];
 end
  save(strcat(dLabels(iSubject).name,'_sigChannel.mat'),'pvalsMCleanPerc');
  %%
   if(sum(pvalsMCleanPerc)>0)
   ieegCarSig = ieegSplit(speechProdChan,:,:);
   ieegCarSig = carFilter(ieegCarSig);
   ieegCarSigBase = ieegCarBase(speechProdChan,:);
   [~,ieegGammaBasedown] = EcogExtractHighGammaTrial(double(squeeze(ieegCarSigBase)),fsD,fsD,[70 150],[0 length(ieegCarSigBase)-1]./fsD,[0 length(ieegCarSigBase)-1]./fsD,[]);
   timeBaseDown = [0:size(ieegGammaBasedown,2)-1]/fsD;
   normFactor = [mean(ieegGammaBasedown(:,timeBaseDown>=5&timeBaseDown<=35),2) std(ieegGammaBasedown(:,timeBaseDown>=5&timeBaseDown<=35),0,2)];
   ieegGamma = [];
    for iTrial = 1:size(ieegCarSig,2)
        iTrial
        [~,ieegGammaTmp] = EcogExtractHighGammaTrial(double(squeeze(ieegCarSig(:,iTrial,:))),fsD,fsD,[70 150],tw,etw,[],[]);
        ieegGamma(:,iTrial,:) = ieegGammaTmp;
    end
    timeGamma = linspace(etw(1),etw(2),size(ieegGamma,3));
    normFactorIntTrial = []
    for iChan = 1:size(ieegCarSig,1)
    normFactorIntTrial(iChan,:) = [mean2(squeeze(ieegGamma(iChan,goodTrialsCommon,timeGamma>=1&timeGamma<1.5))) std2(squeeze(ieegGamma(iChan,goodTrialsCommon,timeGamma>=1&timeGamma<1.5)))];
    end
    ieegGammaNorm = [];
    fsDown = 1000;
    for iTrial = 1:size(ieegCarSig,2)
        iTrial
        [~,ieegGammaTmp] = EcogExtractHighGammaTrial(double(squeeze(ieegCarSig(:,iTrial,:))),fsD,fsDown,[70 150],tw,etwG,normFactorIntTrial,1);
        ieegGammaNorm(:,iTrial,:) = ieegGammaTmp;
    end
    timeGamma = linspace(etwG(1),etwG(2),size(ieegGammaNorm,3));
   %%
trialNames = [];

trialInfoResp = trialInfo(respId);
trialClean = trialInfoResp;
vowel = "aeiou";
ts1 = [-0.1 0];
ts2 = [0 0.1];
timeSelect1 = timeGamma>=ts1(1)&timeGamma<=ts1(2);
timeSelect2 = timeGamma>=ts2(1)&timeGamma<=ts2(2);
phonPairAll = [];
gammaPowerAll = [];
binClass = [];
phonConsClass = [];
phonCVClass = [];
for iTrial = 1:length(trialClean)
    phonPairParse = [];
trialNames{iTrial} = (trialClean(iTrial).sound(1:end-4));
if(strlength(trialNames{iTrial})==3)
    phonPairParse = {trialNames{iTrial}(1:2); trialNames{iTrial}(2:3)};
else
    if(contains(vowel,trialNames{iTrial}(1))&&contains(vowel,trialNames{iTrial}(2)))
    phonPairParse = {trialNames{iTrial}(1:3);trialNames{iTrial}(3:4)};
    else
%         if(contains(vowel,trialNames{iTrial}(1)))
%             phonPairParse = {trialNames{iTrial}(1:2);trialNames{iTrial}(2:4)};
%         end
    phonPairParse = {trialNames{iTrial}(1:2);trialNames{iTrial}(2:4)};
    end
    
end
if(contains(vowel,trialNames{iTrial}(1)))
    binClass = [binClass 1];
    phonConsClass = [phonConsClass 1];
    
%     if(contains("ae",trialNames{iTrial}(1:2)))
%         phonCVClass = [phonCVClass 1];
%     end
    if(contains("a",trialNames{iTrial}(1)))
        phonCVClass = [phonCVClass 1];
    end
    if(contains("iu",trialNames{iTrial}(1)))
        phonCVClass = [phonCVClass 2];
    end
%     if(contains("u",trialNames{iTrial}(1)))
%         phonCVClass = [phonCVClass 3];
%     end
        
    
else
    binClass = [binClass 2];
    if(contains("bpv",trialNames{iTrial}(1)))
        phonConsClass = [phonConsClass 2];
        phonCVClass = [phonCVClass 3];
    end
    if(contains("gk",trialNames{iTrial}(1)))
        phonConsClass = [phonConsClass 3];
        phonCVClass = [phonCVClass 4];
    end
%     if(contains("v",trialNames{iTrial}(1)))
%         phonConsClass = [phonConsClass 2];
%         phonCVClass = [phonCVClass 3];
%     end
end
    
phonPairAll = [phonPairAll; phonPairParse];
gammaPowerPair = [squeeze(mean(ieegGamma(:,iTrial,timeSelect1),3)) squeeze(mean(ieegGamma(:,iTrial,timeSelect2),3))];
gammaPowerAll = [gammaPowerAll gammaPowerPair];

end 
            
%         %%
%       
%         accChan = [];
% %sigChannel = find(pvalsMCleanProdResp);
% matSize = size(ieegGammaNorm);
% gammaFlat = reshape(permute(ieegGammaNorm,[2 1 3]),[matSize(2) matSize(1)*matSize(3)]);
% for iChan = 1:matSize(1)
%     iChan
% gammaFlatBin = [squeeze(ieegGammaNorm(iChan,:,:))];
% phonClassBin = [phonCVClass   ];
% CmatAll = zeros(4,4);
% cvp = cvpartition(phonClassBin,'KFold',5,'Stratify',true);
% accAll = 0;
% labels = {'Low','High','Labial','Dorsal'};
% pComb = nchoosek([1:length(labels)],2);
% CmatCombAll = zeros(size(pComb,1),2,2);
% for nCv = 1:cvp.NumTestSets
%         
%     train = cvp.training(nCv);
%     test = cvp.test(nCv);
%     gTrain = gammaFlatBin(train,:);
%     gTest = gammaFlatBin(test,:);
%     pTrain = phonClassBin(train);
%     pTest = phonClassBin(test);
%     [lossVect,aucVect] = scoreSelect(gTrain,pTrain,50,1,5); % Hyper parameter tuning
% %     figure;
% %     plot(mean(lossVect,1))
%     [~,nDim] = min(mean(lossVect,1)); % Selecting the optimal principal components
%    % mean(squeeze(aucVect(:,nDim,:)),1)
%     [lossMod,Cmat] = pcaDecode(gTrain,gTest,pTrain,...
%                    pTest,nDim);
% accAll = accAll + 1 - lossMod;
% 
% CmatAll = CmatAll + Cmat;
% end
% acc = trace(CmatAll)/sum(sum(CmatAll))
% tp = diag(CmatAll./sum(CmatAll,1));
% [~,maxChan] = max(tp);
% chanLab(iChan) = acc;
% accChan(iChan,:) = tp;
% end

%% 4 - label classification - Low, High, Labial, Dorsal
            %sigChannel = find(pvalsMCleanProdResp);
            matSize = size(ieegGammaNorm);
            gammaFlat = reshape(permute(ieegGammaNorm,[2 1 3]),[matSize(2) matSize(1)*matSize(3)]);
            gammaFlatBin = gammaFlat;
            phonClassBin = phonCVClass;
            CmatAll = zeros(4,4);
            cvp = cvpartition(phonClassBin,'KFold',20,'Stratify',true);
            accAll = 0;
            labels = {'Low','High','Labial','Dorsal'};
            pComb = nchoosek([1:length(labels)],2);
            CmatCombAll = zeros(size(pComb,1),2,2);
            for nCv = 1:cvp.NumTestSets

                train = cvp.training(nCv);
                test = cvp.test(nCv);
                gTrain = gammaFlatBin(train,:);
                gTest = gammaFlatBin(test,:);
                pTrain = phonClassBin(train);
                pTest = phonClassBin(test);
                [lossVect,aucVect] = scoreSelect(gTrain,pTrain,150,1,20); % Hyper parameter tuning
                 figure;
                 plot(mean(lossVect,1))
                [~,nDim] = min(mean(lossVect,1)); % Selecting the optimal principal components
        %        mean(squeeze(aucVect(:,nDim,:)),1)
                [lossMod,Cmat,yhat] = pcaDecode(gTrain,gTest,pTrain,...
                               pTest,nDim);

            accAll = accAll + 1 - lossMod;

            CmatAll = CmatAll + Cmat;
            end
            acc = trace(CmatAll)/sum(sum(CmatAll))
        labels = {'Low','High','Labial','Dorsal'};

        figure; 
        cm = confusionchart(CmatAll,labels,'Normalization','row-normalized');
        sortClasses(cm,labels);
        set(gca,'FontSize',20);
        title(['Accuracy : ' num2str(acc*100) '%']);
   end
   
end
% %% 
% sigPowerAll = {};
% %%
% % specPower = specPower(logical(gridIndex)');
% % pvalsMCleanPerc = pvalsMCleanPerc(logical(gridIndex));
% sigPowerAll{4} = specPower(pvalsMCleanPerc);
% %%
specPowerGrid = specPower(logical(gridIndex));
[~, pvalsMCleanPerc] = fdr(pPerc,0.05);
sigChannelGrid = pvalsMCleanPerc(logical(gridIndex));
sigPowerGrid = specPower(sigChannelGrid);
% %%
% sigPowerStereo = [];
% 
% for iChan = 1:length(sigPowerAll)
%     sigPowerStereo = [sigPowerStereo; sigPowerAll{iChan}];
% end
% %%
% 
% 
figure; distributionPlot({20.*log10(sigPowerStereo),20.*log10(sigPowerGrid), 20.*log10(sigPowerMicro)},'colormap',copper,'showMM',6,'xNames',{'SEEG','macro-grid','\muECoG'},'yLabel','High Gamma Power (dB)');
sigstar({[2,3],[1,3]},[1e-9,1e-31]);