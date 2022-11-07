addpath(genpath('E:\Box Sync\ECoG_Recon\matlab_code\'));
global DUKEDIR
DUKEDIR = 'E:\Box Sync\CoganLab\D_Data\Phoneme_Sequencing';
dLabels = dir(DUKEDIR);
dLabels = dLabels(3:end);

tw = [-3 2]; % time window
etw = [-2.5 1.5];
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
    ieegCAR = []; ieegChan = [];
    ieegCARAll = []; ieegChanAll = [];
    channels = [];
     specRaw = [];
     ieegSplit = [];
     ieegBase = [];
    [ieegSplit,~,trigOnset]=trialIEEGUpdate(Trials,1:length(allChannels),'ResponseStart',tw.*1000);
    ieegSplit = permute(ieegSplit,[2,1,3]);
    %ieegBase = squeeze(trialIEEGUpdate(Trials(1),1:length(allChannels),'ResponseStart',[-40 0].*1000));
    
    
    
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
  %%  
    respId = find(~isnan(trigOnset));
    ieegSplitResp = ieegSplit(:,respId,:);
    trigOns = trigOnset(respId);
    ieegCarImp = carFilterImpedance(ieegSplitResp,badChannels);
    
    
%     ieegCarBase = carFilterImpedance(ieegBase,higRmsId);
%     ieegCarBase = filt60(double(ieegCarBase),fsD);
    %ieegCarImp = carFilterImpedance(ieegSplit,higRmsId);
    ieegCarImpFilt = [];
    for iChan = 1:size(ieegCarImp,1)
        ieegCarImpFilt(iChan,:,:) = filt60(squeeze(double(ieegCarImp(iChan,:,:))),fsD);
    end
    
   % load(strcat(dLabels(iSubject).name,'_sigChannel.mat'))
    time = linspace(tw(1),tw(2),size(ieegCarImpFilt,3));
    
    
   %load(strcat(dLabels(iSubject).name,'_sigChannel.mat'))
    [specRaw,pPerc,pProd]= getSpectrograms(ieegCarImpFilt,[],tw,[-2.5 1.5],[1 200],[-1.5 -1],[-0.25 0.25],[1 3],[50 80],fsD,1);
    [~, pvalsMCleanPerc] = fdr(pPerc,0.05);
    [ir,ic] = numSubplots(length(specRaw));
    
    figure;
    specPower = [];
for iChan = 1:length(specRaw)
     subplot(ir(1),ir(2),iChan);
% specChanMap(spec,[],1:length(channels),etw,[0.7 1.5],iChan)
 %figure;
    [~,specPowerTemp] = specChanMap(specRaw,[],1:length(allChannels),find(pvalsMCleanPerc),[-2.5 1.5],[-1.5 -1],[-0.25 0.25],[1 200],[70 150],[0.7 1.4],iChan);
colormap(parula(4096));
title(allChannels(iChan));
specPower = [specPower; specPowerTemp];
 end
  %save(strcat(dLabels(iSubject).name,'_sigChannel.mat'),'pvalsMCleanPerc');
  %%
  if(sum(pvalsMCleanPerc)>0)
   
   
   ieegCarSig = ieegCarImpFilt(pvalsMCleanPerc,:,:);
   %ieegCarSig = carFilter(ieegCarSig);
   [~,goodtrials] = remove_bad_trials(ieegCarSig,8);
   goodTrialsCommon = extractCommonTrials(goodtrials);
   
  % ieegCarSigBase = ieegCarBase(speechProdChan,:);
  % [~,ieegGammaBasedown] = EcogExtractHighGammaTrial(double(squeeze(ieegCarSigBase)),fsD,fsD,[50 85],[0 length(ieegCarSigBase)-1]./fsD,[0 length(ieegCarSigBase)-1]./fsD,[]);
   %timeBaseDown = [0:size(ieegGammaBasedown,2)-1]/fsD;
  % normFactor = [mean(ieegGammaBasedown(:,timeBaseDown>=5&timeBaseDown<=35),2) std(ieegGammaBasedown(:,timeBaseDown>=5&timeBaseDown<=35),0,2)];
   ieegGamma = [];
    for iTrial = 1:size(ieegCarSig,2)
        iTrial
        [~,ieegGammaTmp] = EcogExtractHighGammaTrial(double(squeeze(ieegCarSig(:,iTrial,:))),fsD,fsD,[70 150],tw,etw,[],[]);
        ieegGamma(:,iTrial,:) = ieegGammaTmp;
    end
    timeGamma = linspace(etw(1),etw(2),size(ieegGamma,3));
    normFactorIntTrial = []
    for iChan = 1:size(ieegCarSig,1)
    normFactorIntTrial(iChan,:) = [mean2(squeeze(ieegGamma(iChan,:,timeGamma>=-1.5&timeGamma<-1))) std2(squeeze(ieegGamma(iChan,:,timeGamma>=-1.5&timeGamma<-1)))];
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
            
        %%
      
        accChan = [];
%sigChannel = find(pvalsMCleanProdResp);
matSize = size(ieegGammaNorm);
gammaFlat = reshape(permute(ieegGammaNorm,[2 1 3]),[matSize(2) matSize(1)*matSize(3)]);
CMatCat = [];
for iChan = 1:matSize(1)
    iChan
gammaFlatBin = [squeeze(ieegGammaNorm(iChan,goodtrials{iChan},timeGamma>=-0.5&timeGamma<=0.5))];
phonClassBin = [phonCVClass(goodtrials{iChan})   ];
CmatAll = zeros(4,4);
cvp = cvpartition(phonClassBin,'KFold',5,'Stratify',true);
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
    [lossVect,aucVect] = scoreSelect(gTrain,pTrain,100,1,5); % Hyper parameter tuning
%     figure;
%     plot(mean(lossVect,1))
    [~,nDim] = min(mean(lossVect,1)); % Selecting the optimal principal components
   % mean(squeeze(aucVect(:,nDim,:)),1)
    [lossMod,Cmat] = pcaDecode(gTrain,gTest,pTrain,...
                   pTest,nDim);
accAll = accAll + 1 - lossMod;

CmatAll = CmatAll + Cmat;
end
acc = trace(CmatAll)/sum(sum(CmatAll))
tp = diag(CmatAll./sum(CmatAll,1));
[~,maxChan] = max(tp);
chanLab(iChan) = acc;
accChan(iChan,:) = tp;
CMatCat(iChan,:,:) = CmatAll;
end
chanLabel = allChannels(pvalsMCleanPerc);
save(strcat(dLabels(iSubject).name,'_indChan_decode.mat'),'chanLab','labels','CMatCat','chanLabel');
% 4 - label classification - Low, High, Labial, Dorsal
            %sigChannel = find(pvalsMCleanProdResp);
            ieegModel = ieegGammaNorm(:,goodTrialsCommon,:);
            matSize = size(ieegModel);
            gammaFlat = reshape(permute(ieegModel,[2 1 3]),[matSize(2) matSize(1)*matSize(3)]);
            %gammaFlatBin = gammaFlat(goodTrialsCommon,:);
            phonClassBin = phonCVClass(goodTrialsCommon);
            CmatAll = zeros(4,4);
            cvp = cvpartition(phonClassBin,'KFold',10,'Stratify',true);
            accAll = 0;
            labels = {'Low','High','Labial','Dorsal'};
            pComb = nchoosek([1:length(labels)],2);
            CmatCombAll = zeros(size(pComb,1),2,2);
            aucAll = zeros(1,length(labels));
            for nCv = 1:cvp.NumTestSets

                train = cvp.training(nCv);
                test = cvp.test(nCv);
                ieegTrain = ieegModel(:,train,timeGamma>=-0.5&timeGamma<=0.5);
                ieegTest = ieegModel(:,test,timeGamma>=-0.5&timeGamma<=0.5);
                matTrain = size(ieegTrain);
                gTrain = reshape(permute(ieegTrain,[2 1 3]),[matTrain(2) matTrain(1)*matTrain(3)]);
                matTest = size(ieegTest);
                gTest = reshape(permute(ieegTest,[2 1 3]),[matTest(2) matTest(1)*matTest(3)]);
                %gTrain = gammaFlatBin(train,:);
                %gTest = gammaFlatBin(test,:);
                pTrain = phonClassBin(train);
                pTest = phonClassBin(test);
                [lossVect,aucVect] = scoreSelect(gTrain,pTrain,100,1,10); % Hyper parameter tuning
                 figure;
                 plot(mean(lossVect,1))
                [~,nDim] = min(mean(lossVect,1)); % Selecting the optimal principal components
        %        mean(squeeze(aucVect(:,nDim,:)),1)
                [lossMod,Cmat,yhat,aucVect] = pcaDecode(gTrain,gTest,pTrain,...
                               pTest,nDim);

            accAll = accAll + 1 - lossMod;
            aucAll = aucAll + aucVect;
            CmatAll = CmatAll + Cmat;
            end
            acc = trace(CmatAll)/sum(sum(CmatAll))
        labels = {'Low','High','Labial','Dorsal'};

        figure; 
        cm = confusionchart(CmatAll,labels,'Normalization','row-normalized');
        sortClasses(cm,labels);
        set(gca,'FontSize',20);
        title(['Accuracy : ' num2str(acc*100) '%']);
        save(strcat(dLabels(iSubject).name,'_decode.mat'),'CmatAll','labels','acc');
   end
end