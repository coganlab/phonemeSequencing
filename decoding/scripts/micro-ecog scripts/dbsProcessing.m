addpath(genpath('E:\Box Sync\ECoG_Recon\matlab_code\'));
global DUKEDIR
DUKEDIR = 'E:\microECoG surgeries\S14 surgery';
dLabels = dir(DUKEDIR);
dLabels = dLabels(3:end);

tw = [-3 2]; % time window
etw = [-2.5 1.5];
etwG = [-1 1];% epoch time window
prtw = [-1.5 -1]; % preonset time window
pstw = [0.25 0.75]; % postonset time window
gammaF = [70 150]; % frequency in Hz
load('channelMap.mat');
selectedChannels = sort(chanMap(~isnan(chanMap)))';
%% data Loading
Experiment = loadExperiment('S14-DBS');
fsD = Experiment.processing.ieeg.sample_rate;
Trials = dbTrials('S14-DBS',Experiment.recording.recording_day,'Speech_OvertMimeMove');

trialFiles = strcat('\',Experiment.recording.recording_day,'\mat\trialInfo.mat');
[ieegSplit,~,trigOnset]=trialIEEGUpdate(Trials,selectedChannels,'ResponseOnset',tw.*1000);
[~,~,trigOffset]=trialIEEGUpdate(Trials,selectedChannels,'ResponseOffset',tw.*1000);
[~,~,trigAuditory]=trialIEEGUpdate(Trials,selectedChannels,'Auditory',tw.*1000);
ieegSplit = permute(ieegSplit,[2,1,3]);

ieegBase = squeeze(trialIEEGUpdate(Trials(end),selectedChannels,'ResponseOnset',[2 42].*1000));
%% Common average referencing
respId = find(~isnan(trigOnset));
ieegSplit = ieegSplit(:,respId,:);
trigOns = trigOnset(respId);
trigOfs = trigOffset(respId);
%%
higImpId = log10(impedance1)>6;
ieegCarImp = carFilterImpedance(ieegSplit,higImpId);
ieegCarBase = carFilterImpedance(ieegBase,higImpId);
ieegCarBase = filt60(double(ieegCarBase),fsD);
ieegCarImpFilt = [];
for iChan = 1:size(ieegCarImp,1)
    ieegCarImpFilt(iChan,:,:) = filt60(squeeze(double(ieegCarImp(iChan,:,:))),fsD);
end
ieegCarImp = ieegCarImpFilt;
%%
[~,goodtrials] = remove_bad_trials(ieegCarImp,10);
[specRawResp,pProd1,pProd2]= getSpectrograms(ieegCarImp,goodtrials,[-3 2],[-2.5 1.5],[1 200],[1 1.5],[-0.9 -0.4],[0.1 0.6],[70 150],fsD,1);
[~, pvalsMCleanProdResp1] = fdr(pProd1,0.05);
[~, pvalsMCleanProdResp2] = fdr(pProd2,0.05);
sigChannel1 = find(pvalsMCleanProdResp1);
sigChannel2 = find(pvalsMCleanProdResp2);
figure;
[~,sigPowerMicro] = specChanMap(specRawResp,chanMap,selectedChannels,find(pvalsMCleanProdResp),[-2.5 1.5],[1 1.5],[-0.25 0.25],[1 200],[70 150],[0.7 1.4],0);

utterOnly = setdiff(sigChannel2,sigChannel1);
figure; plot(timeGamma,squeeze(mean(ieegGamma(utterOnly,:,:),2)),'color',[0 0 0]+0.5);
figure; plot(timeGamma,squeeze(mean(ieegGamma(sigChannel2,:,:),2)),'color',[0 0 0]+0.5);

figure; plot(timeGamma,mean(squeeze(mean(ieegGamma(utterOnly,:,:),2)),1));
hold on;
plot(timeGamma,mean(squeeze(mean(ieegGamma(sigChannel2,:,:),2)),1));
%%
%%
load('sigChannel.mat');
fsDown =100;
[~,goodtrials] = remove_bad_trials(ieegCarImp,12);
goodTrialsCommon = extractCommonTrials(goodtrials(pvalsMCleanProdResp));
ieegCarClean = ieegCarImp(:,goodTrialsCommon,:);

[~,ieegGammaBasedown] = EcogExtractHighGammaTrial(double(squeeze(ieegCarBase)),fsD,fsDown,[70 150],[0 length(ieegCarBase)-1]./fsD,[0 length(ieegCarBase)-1]./fsD,[]);
   timeBaseDown = [0:size(ieegGammaBasedown,2)-1]/fsDown;
   normFactor = [mean(ieegGammaBasedown(:,timeBaseDown>=5&timeBaseDown<=35),2) std(ieegGammaBasedown(:,timeBaseDown>=5&timeBaseDown<=35),0,2)];
  
ieegGamma = [];
for iTrial = 1:size(ieegCarClean,2)
    iTrial
    [~,ieegGammaTmp] = EcogExtractHighGammaTrial(double(squeeze(ieegCarClean(:,iTrial,:))),fsD,fsDown,[70 150],tw,etwG,normFactor,2);
    ieegGamma(:,iTrial,:) = ieegGammaTmp;
end
timeGamma = linspace(etwG(1),etwG(2),size(ieegGamma,3));
ieegGammaClean = ieegGamma;

%%

    allGammaUtter = squeeze(mean(ieegGamma(utterOnly,:,:),2));
    meanGammaUtter = mean(allGammaUtter,1);
    allGammaPlan = squeeze(mean(ieegGamma(sigChannel2,:,:),2));
    meanGammaPlan = mean(allGammaPlan,1);
    figure;
    plot(timeGamma,allGammaUtter','color',[0 0 1 0.2]);
    hold on;
    plot(timeGamma,meanGammaUtter,'color',[0 0 1]);
    plot(timeGamma,allGammaPlan','color',[1 0 0 0.2]);
    hold on;
    plot(timeGamma,meanGammaPlan,'color',[1 0 0]);
%%
trialNames = [];

trialInfoResp = trialInfo(respId);
trialClean = trialInfoResp(goodTrialsCommon);
vowel = "aeiou";
ts1 = [-1 -0.5];
ts2 = [0 0.5];
timeSelect1 = timeGamma>=ts1(1)&timeGamma<=ts1(2);
timeSelect2 = timeGamma>=ts2(1)&timeGamma<=ts2(2);
phonPairAll = [];
gammaPowerAll = [];
binClass = [];
phonConsClass = [];
phonCVClass = [];
phonCvClassEnd = [];
phonIndClass = [];
phonIndClassEnd = [];
for iTrial = 1:length(trialClean)
    phonPairParse = [];
trialNames{iTrial} = (trialClean{iTrial}.sound(1:end-4));
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
        if(contains("e",trialNames{iTrial}(2)))
        phonIndClass = [phonIndClass 2];
        else
            phonIndClass = [phonIndClass 1];
        end
    end
    if(contains("iu",trialNames{iTrial}(1)))
        phonCVClass = [phonCVClass 2];
        if(contains("i",trialNames{iTrial}(1)))
            phonIndClass = [phonIndClass 3];
        end
        if(contains("u",trialNames{iTrial}(1)))
            phonIndClass = [phonIndClass 4];
        end
    end
    
     if(contains("ae",trialNames{iTrial}(end)))
        phonCvClassEnd = [phonCvClassEnd 1];
        if(contains("e",trialNames{iTrial}(end)))
        phonIndClassEnd = [phonIndClassEnd 2];
        else
            phonIndClassEnd = [phonIndClassEnd 1];
        end
    end
    if(contains("iu",trialNames{iTrial}(end)))
        phonCvClassEnd = [phonCvClassEnd 2];
        if(contains("i",trialNames{iTrial}(end)))
            phonIndClassEnd = [phonIndClassEnd 3];
        end
        if(contains("u",trialNames{iTrial}(end)))
            phonIndClassEnd = [phonIndClassEnd 4];
        end
    end
%     if(contains("u",trialNames{iTrial}(1)))
%         phonCVClass = [phonCVClass 3];
%     end
        
    
else
    binClass = [binClass 2];
    if(contains("bpv",trialNames{iTrial}(1)))
        phonConsClass = [phonConsClass 2];
        phonCVClass = [phonCVClass 3];
        if(contains("b",trialNames{iTrial}(1)))
            phonIndClass = [phonIndClass 5];
        end
        if(contains("p",trialNames{iTrial}(1)))
            phonIndClass = [phonIndClass 6];
        end
        if(contains("v",trialNames{iTrial}(1)))
            phonIndClass = [phonIndClass 7];
        end
    end
    if(contains("gk",trialNames{iTrial}(1)))
        phonConsClass = [phonConsClass 3];
        phonCVClass = [phonCVClass 4];
        if(contains("g",trialNames{iTrial}(1)))
            phonIndClass = [phonIndClass 8];
        end
        if(contains("k",trialNames{iTrial}(1)))
            phonIndClass = [phonIndClass 9];
        end
    end
    
    if(contains("bpv",trialNames{iTrial}(end)))
        %phonConsClass = [phonConsClass 2];
        phonCvClassEnd = [phonCvClassEnd 3];
        
        if(contains("b",trialNames{iTrial}(end)))
            phonIndClassEnd = [phonIndClassEnd 5];
        end
        if(contains("p",trialNames{iTrial}(end)))
            phonIndClassEnd = [phonIndClassEnd 6];
        end
        if(contains("v",trialNames{iTrial}(end)))
            phonIndClassEnd = [phonIndClassEnd 7];
        end
    end
    if(contains("gk",trialNames{iTrial}(end)))
        %phonConsClass = [phonConsClass 3];
        phonCvClassEnd = [phonCvClassEnd 4];
        if(contains("g",trialNames{iTrial}(end)))
            phonIndClassEnd = [phonIndClassEnd 8];
        end
        if(contains("k",trialNames{iTrial}(end)))
            phonIndClassEnd = [phonIndClassEnd 9];
        end
    end
%     
end
    
phonPairAll = [phonPairAll; phonPairParse];
gammaPowerPair = [squeeze((ieegGamma(:,iTrial,timeSelect1))) ;squeeze((ieegGamma(:,iTrial,timeSelect2)))];
gammaPowerAll = [gammaPowerAll gammaPowerPair];

end
gammaPowerSeries = [ieegGamma(:,:,timeSelect1) ieegGamma(:,:,timeSelect2)];
%% first phoneme classification
sigChannel = find(pvalsMCleanProdResp);
matSize = size(ieegGammaClean(sigChannel,:,:));
ieegGammaChan = ieegGammaClean(sigChannel,:,timeGamma>=-0.75&timeGamma<=0.75);

gammaFlat1 = reshape(permute(ieegGammaClean(sigChannel,:,:),[2 1 3]),[matSize(2) matSize(1)*matSize(3)]);
gammaFlat2 = reshape(permute(ieegGammaClean(sigChannel,:,:),[2 1 3]),[matSize(2) matSize(1)*matSize(3)]);
for iSelect = 1
[pTrainPer,perId] = datasample(phonIndClass,round(1*length(phonIndClass)),'Replace',false);
gTrainPer = gammaFlat1(perId,:);
vowelid = find(pTrainPer==1 | pTrainPer==2 | pTrainPer==3 | pTrainPer==4);
consid = find(pTrainPer==5 | pTrainPer==6 | pTrainPer==7 | pTrainPer==8 | pTrainPer==9);
gammaFlat = [gTrainPer];
binClass = [pTrainPer];
phonUnique = unique(binClass);
CmatAll = zeros(length(phonUnique),length(phonUnique));
%nFold = 3;
%cvp = cvpartition(binClass,'KFold',nFold,'Stratify',true);
cvp = cvpartition(binClass,'LeaveOut');
accAll = 0;
ytest = []; ypred = [];
[lossVect,aucVect] = scoreSelect(gammaFlat,binClass,42,0,10); % Hyper parameter tuning
    [~,nDim] = min(mean(lossVect,1));
for nCv = 1:cvp.NumTestSets
        
    train = cvp.training(nCv);
    test = cvp.test(nCv);
    ytest = [ytest binClass(test)];
%     [lossVect,aucVect] = scoreSelect(gammaFlatSeries(train,:),binClass(train),30); % Hyper parameter tuning
%     [~,nDim] = min(mean(lossVect,1));
    [lossMod,Cmat,ypredTemp] = pcaDecode(gammaFlat(train,:),gammaFlat(test,:),binClass(train),...
                   binClass(test),nDim);
               accAll = accAll + 1 - lossMod;
    ypred = [ypred ypredTemp'];
%CmatAll = CmatAll + Cmat;
end
CmatAll = confusionmat(ytest,ypred);
labels = {'a','ae','i','u','b','p','v','g','k'};
%labels = {'a','ae','i','u'};
%labels = {'b','p','v','g','k'};
acc = trace(CmatAll)/sum(sum(CmatAll))
figure; 
cm = confusionchart(CmatAll,labels,'Normalization','row-normalized');
sortClasses(cm,labels);
set(gca,'FontSize',15);
title(['Accuracy :' num2str(acc*100)]);
end


CmatAll = confusionmat(ytestall,ypredall);
labels = {'a','ae','i','u','b','p','v','g','k'};
%labels = {'a','ae','i','u'};
%labels = {'b','p','v','g','k'};
acc = trace(CmatAll)/sum(sum(CmatAll))
figure; 
cm = confusionchart(CmatAll,labels,'Normalization','row-normalized');
sortClasses(cm,labels);
set(gca,'FontSize',15);
title(['Accuracy :' num2str(acc*100)]);
%% Binary classifier
sigChannel = find(pvalsMCleanProdResp);
matSize = size(ieegGammaClean(sigChannel,:,:));
gammaFlat = reshape(permute(ieegGammaClean(sigChannel,:,:),[2 1 3]),[matSize(2) matSize(1)*matSize(3)]);

CmatAll = zeros(2,2);
cvp = cvpartition(binClass,'KFold',5,'Stratify',true);
accAll = 0;
for nCv = 1:cvp.NumTestSets
        
    train = cvp.training(nCv);
    test = cvp.test(nCv);
    [lossVect,aucVect] = scoreSelect(gammaFlat(train,:),binClass(train),30); % Hyper parameter tuning
    figure;
    plot(mean(lossVect,1))
    [~,nDim] = min(mean(lossVect,1)); % Selecting the optimal principal components
    mean(squeeze(aucVect(:,nDim,:)),1)
    [lossMod,Cmat] = pcaDecode(gammaFlat(train,:),gammaFlat(test,:),binClass(train),...
                   binClass(test),nDim);

accAll = accAll + 1 - lossMod;

CmatAll = CmatAll + Cmat;
end

acc = trace(CmatAll)/sum(sum(CmatAll))

labels = {'Vowels','Consonants'};
figure; 
cm = confusionchart(CmatAll,labels);
sortClasses(cm,labels);
set(gca,'FontSize',15);
title(['Accuracy :' num2str(acc*100)]);

%% Phoneme Consonant classifier
sigChannel = find(pvalsMCleanProdResp);
matSize = size(ieegGammaClean(sigChannel,:,:));
gammaFlat = reshape(permute(ieegGammaClean(sigChannel,:,:),[2 1 3]),[matSize(2) matSize(1)*matSize(3)]);
CmatAll = zeros(3,3);
cvp = cvpartition(phonConsClass,'KFold',5,'Stratify',true);
accAll = 0;
for nCv = 1:cvp.NumTestSets
        
    train = cvp.training(nCv);
    test = cvp.test(nCv);
    [lossVect,aucVect] = scoreSelect(gammaFlat(train,:),phonConsClass(train),22); % Hyper parameter tuning
    figure;
    plot(mean(lossVect,1))
    [~,nDim] = min(mean(lossVect,1)); % Selecting the optimal principal components
    mean(squeeze(aucVect(:,nDim,:)),1)
    [lossMod,Cmat] = pcaDecode(gammaFlat(train,:),gammaFlat(test,:),phonConsClass(train),...
                   phonConsClass(test),nDim);

accAll = accAll + 1 - lossMod;

CmatAll = CmatAll + Cmat;
end

acc = trace(CmatAll)/sum(sum(CmatAll))

labels = {'Vowels','Labial','Dorsal Tongue'};

figure; 
cm = confusionchart(CmatAll,labels);
sortClasses(cm,labels);
set(gca,'FontSize',15);
title(['Accuracy :' num2str(acc*100)]);
%% Phoneme Consonant Vowel classifier
sigChannel = find(pvalsMCleanProdResp);
matSize = size(ieegGammaClean(sigChannel,:,:));
gammaFlat = reshape(permute(ieegGammaClean(sigChannel,:,:),[2 1 3]),[matSize(2) matSize(1)*matSize(3)]);
% matSizeParse = size(gammaPowerSeries(sigChannel,:,:));
% gammaFlatSeries = reshape(permute(gammaPowerSeries(sigChannel,:,:),[2 1 3]),[matSizeParse(2) matSizeParse(1)*matSizeParse(3)]);
gammaFlatBin = [ gammaFlat];
phonClassBin = [ phonCVClass ];
CmatAll = zeros(4,4);
cvp = cvpartition(phonClassBin,'KFold',15,'stratify',true);
accAll = 0;
labels = {'Low','High','Labial','Dorsal'};
pComb = nchoosek([1:length(labels)],2);
CmatCombAll = zeros(size(pComb,1),2,2);
accAllPer = [];
for nCv = 1:cvp.NumTestSets
        
    train = cvp.training(nCv);
    test = cvp.test(nCv);
    gTrain = gammaFlatBin(train,:);
    gTest = gammaFlatBin(test,:);
    pTrain = phonClassBin(:,train);
    pTest = phonClassBin(:,test);
    
%     figure;
%     plot(mean(lossVect,1))
accPer = [];
for trainPer = 1:4
    [pTrainPer,perId] = datasample(pTrain,round(trainPer/4*length(pTrain)),'Replace',false);
    gTrainPer = gTrain(perId,:);
    [lossVect,aucVect] = scoreSelect(gTrainPer,pTrainPer,20); % Hyper parameter tuning
    [~,nDim] = min(mean(lossVect,1)); % Selecting the optimal principal components
    %mean(squeeze(aucVect(:,nDim,:)),1)
    [lossMod,Cmat,sTrain,sTest] = pcaDecode(gTrainPer,gTest,pTrainPer,...
                   pTest,nDim);
    accPer = [accPer 1-lossMod];
end
accAllPer = [accAllPer; accPer];
                 
%                
%     
%     for bClassId = 1 : size(pComb,1)
%         
%         pTrainCombId = pTrain==pComb(bClassId,1)|pTrain==pComb(bClassId,2);
%         pTestCombId = pTest==pComb(bClassId,1)|pTest==pComb(bClassId,2);
%         sTrainComb = sTrain(pTrainCombId,:);
%         sTestComb = sTest(pTestCombId,:);
%         pTrainComb = pTrain(pTrainCombId);
%         pTestComb = pTest(pTestCombId);
%         
%         linModel = fitcdiscr(sTrainComb,pTrainComb,'DiscrimType','linear','CrossVal','off');      
%         [yhat,yscore] = predict(linModel,sTestComb);        
%         CmatTemp =  confusionmat(pTestComb,yhat);
%         CmatCombAll(bClassId,:,:) = squeeze(CmatCombAll(bClassId,:,:)) + CmatTemp;
%     end
%     
               
    

accAll = accAll + 1 - lossMod;

CmatAll = CmatAll + Cmat;
end
%
acc = trace(CmatAll)/sum(sum(CmatAll))
labels = {'Low','High','Labial','Dorsal'};

figure; 
cm = confusionchart(CmatAll,labels,'Normalization','row-normalized');
sortClasses(cm,labels);
set(gca,'FontSize',20);
title(['Accuracy : ' num2str(acc*100) '%']);
%%
accMat = nan(4,4);
accAll = [];
for bClassId = 1 : size(pComb,1)
    cmatcomb = squeeze(CmatCombAll(bClassId,:,:));
    acc = trace(cmatcomb)/sum(sum(cmatcomb));
    labelsTemp = {labels{pComb(bClassId,1)},labels{pComb(bClassId,2)} };
    accMat(pComb(bClassId,1),pComb(bClassId,2)) = acc;
    accMat(pComb(bClassId,2),pComb(bClassId,1)) = acc;
    accAll = [accAll acc];
%     figure; 
% cm = confusionchart(cmatcomb,labelsTemp,'Normalization','row-normalized');
% sortClasses(cm,labelsTemp);
% set(gca,'FontSize',20);
% title(['Accuracy : ' num2str(acc*100) '%']);
end

figure; 
accIm = accMat.*100;
ii = ones(size(accIm));
idx = triu(ii,1);
accIm(~idx) = nan;
b = imagesc(accIm);

set(b,'AlphaData',~isnan(accIm));
%insertMarker(accMat,[1 4],'color','red')
caxis([50 100])
colormap(flipud(gray(4096)));
colorbar;
set(gca,'YDir', 'normal');
set(gca,'xtick',1:4)
set(gca,'ytick',1:4);

set(gca,'YTickLabel',labels,'YAxisLocation','left');
hold on;
set(gca,'XTickLabel',round(nansum(accMat)*100/3),'XAxisLocation','top');
set(gca,'XTickLabel',labels,'XAxisLocation','bottom');
set(gca,'FontSize',20);
%% Individual Classification
accChan = [];
chanLab = [];
sigChannel = find(pvalsMCleanProdResp);
matSize = size(ieegGammaClean);
gammaFlat = reshape(permute(ieegGammaClean,[2 1 3]),[matSize(2) matSize(1)*matSize(3)]);
ts = -0.75;
for iTime = 1:length(ts)
timeSelect = timeGamma>=ts(iTime)&timeGamma<=ts(iTime)+1.5;
for iChan = 1:matSize(1)
    iChan
gammaFlatBin = [squeeze(ieegGammaClean(iChan,:,timeSelect))];
phonClassBin = [phonCVClass   ];
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
    [lossVect,aucVect] = scoreSelect(gTrain,pTrain,22,1,5); % Hyper parameter tuning
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
chanLab(iTime,iChan) = acc;
accChan(iTime,iChan,:) = tp;
end
end
for iTime = 1:16
figure;
for iPhon = 1:size(accChan,3)
    subplot(2,2,iPhon)
    chanView(squeeze(accChan(iTime,:,iPhon)).*100,chanMap,selectedChannels,labels{iPhon},[25 100]);
    colormap(flipud(bone(4096)));
    set(gca,'FontSize',20);
    set(gca,'xtick',[])
    set(gca,'ytick',[])
   
end
sgtitle(num2str(ts(iTime)));
end
%% Gamma power average
phonUniq = unique(phonCVClass);
timeSelect = timeGamma>=-1&timeGamma<=-0.5;
for iPhon = phonUniq
    for iChan = 1:128
    gammaPowerAverage(iChan,iPhon) = mean2(squeeze(ieegGammaClean(iChan,phonCVClass==iPhon,timeSelect)));
    end
end
for iPhon = phonUniq
    
    gammaPowerTime(:,iPhon,:) = mean((ieegGammaClean(:,phonCVClass==iPhon,:)),2);
   
end
figure;
for iPhon = 1:size(accChan,2)
    subplot(2,2,iPhon)
    chanView(gammaPowerAverage(:,iPhon)',chanMap,selectedChannels,labels{iPhon},[-2 2]);
    colormap((jet(4096)));
    set(gca,'FontSize',20);
    set(gca,'xtick',[])
    set(gca,'ytick',[])
  
end
sigchannel = setdiff(1:128, find(pvalsMCleanProdResp));
figure;
C = {'k','b','r','g'};
for iPhon = 1:4
    meanGamma = squeeze(mean(gammaPowerTime(sigchannel,iPhon,:),1))';
    stdGamma = squeeze(std(gammaPowerTime(sigchannel,iPhon,:),0,1))';

     H = shadedErrorBar(timeGamma,meanGamma,stdGamma,'lineprops',C{iPhon});
     H.mainLine.LineWidth = 2;
     H.edge(1).LineWidth = 1;
     H.edge(2).LineWidth = 1;
     set(gca,'FontSize',20);
        hold on;
        
        
end
figure;
colorPhon = [0 0 1; 0 0 0; 1 0 0; 0 1 0];
for iPhon = 1:4
    meanGamma = squeeze(mean(gammaPowerTime(sigchannel,iPhon,:),1))';
 plot(timeGamma,meanGamma','color',colorPhon(iPhon,:),'LineWidth',2);
 hold on;
       
        set(gca,'FontSize',20);
        end
for iPhon = 1:4
    meanGamma = squeeze(mean(gammaPowerTime(sigchannel,iPhon,:),1))';
    allGamma = squeeze(gammaPowerTime(sigchannel,iPhon,:))';
       ps =  plot(timeGamma,allGamma','color',[colorPhon(iPhon,:) 0.2]);
       hold on;    
end
legend(labels);
for iPhon = 1:4
    meanGamma = squeeze(mean(gammaPowerTime(sigchannel,iPhon,:),1))';
 plot(timeGamma,meanGamma','color',colorPhon(iPhon,:),'LineWidth',2);
 hold on;
       
        set(gca,'FontSize',20);
        end

%% direct discriminant analysis
sigChannel = find(pvalsMCleanProdResp);
matSize = size(ieegGammaClean(sigChannel,:,:));
gammaFlat = reshape(permute(ieegGammaClean(sigChannel,:,:),[2 1 3]),[matSize(2) matSize(1)*matSize(3)]);
matSizeParse = size(gammaPowerSeries(sigChannel,:,:));
gammaFlatSeries = reshape(permute(gammaPowerSeries(sigChannel,:,:),[2 1 3]),[matSizeParse(2) matSizeParse(1)*matSizeParse(3)]);
gammaFlatBin = [ gammaFlat];
phonClassBin = [ phonCVClass ];
CmatAll = zeros(4,4);
cvp = cvpartition(length(phonClassBin),'KFold',4);
accAll = 0;
labels = {'Low','High','Labial','Dorsal'};
pComb = nchoosek([1:length(labels)],2);
CmatCombAll = zeros(size(pComb,1),2,2);
for nCv = 1:cvp.NumTestSets
        
    train = cvp.training(nCv);
    test = cvp.test(nCv);
    gTrain = gammaFlatBin(train,:);
    gTest = gammaFlatBin(test,:);
    meanTrain = mean(gTrain,1);
    stdTrain = std(gTrain,0,1);
    pTrain = phonClassBin(:,train);
    pTest = phonClassBin(:,test);
    gTrainNorm = (gTrain - meanTrain)./stdTrain;
    gTestNorm = (gTest - meanTrain)./stdTrain;
    phonModel = fitcdiscr(gTrainNorm,pTrain);
    acc = 1 - loss(phonModel,gTestNorm,pTest)
end
%% Decision tree feature importance
labels = {'Low','High','Labial','Dorsal'};
load('channelMap.mat')
selectedChannels = 1:128;
figure;
featureImportance = [];
for iFold = 1:length(cvResult)
    featureImportance(iFold,:,:) = cvResult{iFold}.featImp;
    auc(iFold) = mean(cvResult{iFold}.AUC);
end
goodMod = find(auc>0.7);
figure;
for iPhon = 1:4
    subplot(2,2,iPhon)
    chanView( (log(squeeze(mean((featureImportance(goodMod,iPhon,:)),1)))),chanMap,selectedChannels,labels{iPhon},[-5.5 -2.5]);
    colormap(flipud(bone(4096)));
    set(gca,'FontSize',20);
    set(gca,'xtick',[])
    set(gca,'ytick',[])
  
end

%%
load('sigChannel.mat');
phonSeqModel = fitcdiscr(zscore(gammaAllOccur(pvalsMCleanProdResp,:)'),pLabelOccur,'KFold',5);
%cvphonSeqModel = crossval(phonSeqModel);
acc = 1-kfoldLoss(phonSeqModel)
ypred = kfoldPredict(phonSeqModel);
figure;
confusionchart(pLabelOccur,ypred);
 
%% STRF Audio preprocessing
% 
% 
% micAudPath = 'E:\microECoG surgeries\S14 surgery\audioFiles';
% [micClean,fsMic] = microphoneMergeTask(micAudPath,microphone);
micSplit = squeeze(splitIeeg(micCleanAudio',trigOns.*fsMic/fsD,etwG,fsMic));
timeMic = linspace(etwG(1),etwG(2),size(micSplit,2));

%% melSpectrogram 
SMel = [];
for iTrial = 1:size(micSplit,1)
[SMel(iTrial,:,:),cf,tmel] = melSpectrogram(micSplit(iTrial,:)',fsMic,'WindowLength',200,...
                                'OverlapLength',100,'NumBands',64,'FrequencyRange',[2 6e3],...
                                'FFTLength',8192);
                            
%SMel(iTrial,:,:) =  log10(squeeze(SMel(iTrial,:,:)./mean(SMel(iTrial,:,tmel>=2),3)))     ; 
SMel(iTrial,:,:) = log10(squeeze(SMel(iTrial,:,:)));
figure;
                            
imagesc(tmel,cf,(squeeze(SMel(iTrial,:,:))));
%caxis([-15 -4]);
set(gca,'YDir', 'normal');
end

%% Clean signals
ieegGammaClean = ieegGamma;
SMelClean = SMel(goodTrialsCommon,:,:);
%% 
%[audioPathSort,sortId] = sort(audioPathUsed(goodTrialsCommon)');
[train,test] = crossvalind('HoldOut',length(goodTrialsCommon),0.3);
%cvp = cvpartition(length(goodTrialsCommon),'HoldOut',0.3);
% train = sortId(13:48);
% test = sortId(1:12);
% train = train(randperm(length(train)));
% test = test(randperm(length(test)));


XTrainSequence = ieegGammaClean(:,train,:);
YTrainSequence = SMelClean(train,:,:);
XTestSequence = ieegGammaClean(:,test,:);
YTestSequence = SMelClean(test,:,:);

% Spectrogram Normalization
for iFreq = 1:size(YTrainSequence,2)
    YTrSFreq = squeeze(YTrainSequence(:,iFreq,:));
    YTeSFreq = squeeze(YTestSequence(:,iFreq,:));
    meanFreq = mean2(YTrSFreq);
    stdFreq = std2(YTrSFreq);
    YTrNorm = (YTrSFreq - meanFreq)./stdFreq;
    YTeNorm = (YTeSFreq - meanFreq)./stdFreq;
    YTrainSequence(:,iFreq,:) = reshape(YTrNorm,[size(YTrNorm,1),1,size(YTrNorm,2)]);
    YTestSequence(:,iFreq,:) = reshape(YTeNorm,[size(YTeNorm,1),1,size(YTeNorm,2)]);
end

for iChan = 1:size(XTrainSequence,1)
    XTrainChan = squeeze(XTrainSequence(iChan,:,:));
    XTestChan = squeeze(XTestSequence(iChan,:,:));
    meanChan = mean2(XTrainChan);
    stdChan = std2(XTrainChan);
    XTrainNorm = (XTrainChan-meanChan)./stdChan;
    XTestNorm = (XTestChan-meanChan)./stdChan;
    XTrainSequence(iChan,:,:) = XTrainNorm;
    XTestSequence(iChan,:,:) = XTestNorm;
end


XTrainSequence = XTrainSequence(:,:,1:size(YTrainSequence,3));
XTestSequence = XTestSequence(:,:,1:size(YTestSequence,3));
save('EncodeModel1.mat','XTrainSequence','YTrainSequence','XTestSequence','YTestSequence','cf');
%%
for iChan = 1:size(strfCoef,1)
figure;
    tspec = delays;
        strf2visualize =  squeeze(strfCoef(iChan,:,:));
        %subplot(size(chanMap,1),size(chanMap,2),find(ismember(chanMap',selectedChannels(i))));
               imagesc(tspec,[],strf2visualize);
                set(gca,'YDir', 'normal');
                caxis([-max(max(abs(strf2visualize))) max(max(abs(strf2visualize)))]);
        %if(max(max(abs(strf2visualize))) ==0)
        %tvimage((strf2visualize'),'XRange',[tspec(1),tspec(end)],'YRange',[cf(1),cf(end)]);
         set(gca,'ytick',1:10:length(cf))
        set(gca,'YTickLabel',round(cf(1:10:end))./1000)
        xlabel('Time-advance (s)');
        ylabel('Mel - Frequency (kHz)');
        set(gca,'FontSize',15);
        %else
        %tvimage((strf2visualize'),'XRange',[tspec(1),tspec(end)],'YRange',[cf(1),cf(end)],'CLim', [-max(max(abs(strf2visualize))) max(max(abs(strf2visualize)))]);
        %end
%         axis off;
%         set(gca,'xtick',[],'ytick',[])
        colormap(jet(4096));
end
%         figure
%         
%         chanView(mean(spearmanCorrelation,2)',chanMap,selectedChannels,'Prediction score',[0 0.3])