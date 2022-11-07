fsD = 2000;
fsDown = 200;
load('channelMap.mat')
load('sigChannel.mat');
ieegBase = (h5read('S14_maxlike_frames.h5','/baseline_cleaned'))';
ieegSplit = h5read('S14_maxlike_frames.h5','/trials_cleaned');
ieegSplit = permute(ieegSplit,[2,3,1]);
[~,ieegGammaBasedown] = EcogExtractHighGammaTrial(double(squeeze(ieegBase)),fsD,fsDown,[75 150],[0 length(ieegBase)-1]./fsD,[0 length(ieegBase)-1]./fsD,[]);
timeBaseDown = [0:size(ieegGammaBasedown,2)-1]/fsDown;
normFactor = [mean(ieegGammaBasedown(:,timeBaseDown>=10&timeBaseDown<=40),2) std(ieegGammaBasedown(:,timeBaseDown>=10&timeBaseDown<=40),0,2)];

tw = [-3 2]; % time window
etw = [-2.5 1.5];
etwG = [-1 1];
channel = 1:128;
fsDown =200;
[~,goodtrials] = remove_bad_trials(ieegSplit,12);
goodTrialsCommon = extractCommonTrials(goodtrials(pvalsMCleanProdResp));
ieegCarClean = ieegSplit(:,goodTrialsCommon,:);

ieegGamma = [];
for iTrial = 1:size(ieegCarClean,2)
    iTrial
    [~,ieegGammaTmp] = EcogExtractHighGammaTrial(double(squeeze(ieegCarClean(:,iTrial,:))),fsD,fsD,[70 150],tw,etwG,normFactor);
    ieegGamma(:,iTrial,:) = ieegGammaTmp;
end
timeGamma = linspace(etwG(1),etwG(2),size(ieegGamma,3));
ieegGammaClean = ieegGamma;
travelGamma = squeeze(mean(ieegGammaClean,2));
travellingWaveMovie(travelGamma,fsD,chanMap,channel,timeGamma,[-1 1],[-2 2],240,'S14gammaInducedResponse2');
%timeGamma = linspace(tw(1),tw(2),size(ieegGammaNorm,3));

%%
load('responseId.mat');
load('S14_Block_3_TrialData.mat');
trialNames = [];

trialInfoResp = trialInfo(respId);
trialClean = trialInfoResp(goodTrialsCommon);
vowel = "aeiou";
ts1 = [-1 0.25];
ts2 = [-0.25 1];
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

%%
%% first phoneme classification
sigChannel = find(pvalsMCleanProdResp);
matSize = size(ieegGammaClean(:,:,:));
gammaFlat1 = reshape(permute(ieegGammaClean(:,:,:),[2 1 3]),[matSize(2) matSize(1)*matSize(3)]);
gammaFlat2 = reshape(permute(ieegGammaClean(:,:,:),[2 1 3]),[matSize(2) matSize(1)*matSize(3)]);
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
%%
CmatAll = confusionmat(ytest,ypredict);
labels = {'Low','High','Labial','Dorsal'};
%labels = {'a','ae','i','u'};
%labels = {'b','p','v','g','k'};
acc = trace(CmatAll)/sum(sum(CmatAll))
figure; 
cm = confusionchart(CmatAll,labels,'Normalization','row-normalized');
sortClasses(cm,labels);
set(gca,'FontSize',15);
title(['Accuracy :' num2str(acc*100)]);