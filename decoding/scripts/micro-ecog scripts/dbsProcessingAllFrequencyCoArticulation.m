addpath(genpath('E:\Box Sync\ECoG_Recon\matlab_code\'));
global DUKEDIR
DUKEDIR = 'E:\microECoG surgeries\S14 surgery';
dLabels = dir(DUKEDIR);
dLabels = dLabels(3:end);

tw = [-2 2]; % time window
etw = [-1.5 1.5];
etwG = [-1 1];% epoch time window
prtw = [-1.5 -1]; % preonset time window
pstw = [0.25 0.75]; % postonset time window
gammaF = [70 150]; % frequency in Hz
load('channelMap.mat');
selectedChannels = sort(chanMap(~isnan(chanMap)))';
Experiment = loadExperiment('S14');
fsD = Experiment.processing.ieeg.sample_rate;
%%
Trials = dbTrials('S14',Experiment.recording.recording_day,'Speech_OvertMimeMove');

trialFiles = strcat('\',Experiment.recording.recording_day,'\mat\trialInfo.mat');
[ieegSplit1,~,trigOnset]=trialIEEGUpdate(Trials,selectedChannels,'phon1Onset','ieeg',tw.*1000);
ieegSplit1 = permute(ieegSplit1,[2,1,3]);
psdAll = [];
for iChan = 1:size(ieegSplit1,1)
    [psdAll(iChan,:,:),f] = pwelch(squeeze(ieegSplit1(iChan,:,:))',2*fsD,1*fsD,4096,fsD);
end

[ieegSplit2,~,~]=trialIEEGUpdate(Trials,selectedChannels,'phon2Onset','ieeg',tw.*1000);
%[ieegSplit3,~,~]=trialIEEGUpdate(Trials,selectedChannels,'phon3Onset',tw.*1000);
[~,~,trigOffset]=trialIEEGUpdate(Trials,selectedChannels,'ResponseOffset','ieeg',tw.*1000);
[~,~,trigAuditory]=trialIEEGUpdate(Trials,selectedChannels,'Auditory','ieeg',tw.*1000);

ieegSplit2 = permute(ieegSplit2,[2,1,3]);
%ieegSplit3 = permute(ieegSplit3,[2,1,3]);
phon1 = [Trials.phon1];
phon2 = [Trials.phon2];
phon3 = [Trials.phon3];
ieegBase = squeeze(trialIEEGUpdate(Trials(end),selectedChannels,'ResponseOnset','ieeg',[2 42].*1000));
 
respId = find(~isnan(trigOnset));
ieegSplit1 = ieegSplit1(:,respId,:);
ieegSplit2 = ieegSplit2(:,respId,:);
%ieegSplit3 = ieegSplit3(:,respId,:);
phon1 = phon1(respId);
phon2 = phon2(respId);
phon3 = phon3(respId);
phon12 = append(phon1,phon2);
phon23 = append(phon2,phon3);
%%
ieegSplitAll = cat(2,ieegSplit1, ieegSplit2);
phonAll = [phon12 phon23];
trigOns = trigOnset(respId);
trigOfs = trigOffset(respId);

%% Common average referencing
higImpId = log10(impedance1)>6;
ieegCarImp = carFilterImpedance(ieegSplitAll,higImpId);
ieegCarBase = carFilterImpedance(ieegBase,higImpId);
ieegCarBase = filt60(double(ieegCarBase),fsD);
% ieegCarImpFilt = [];
% for iChan = 1:size(ieegCarImp,1)
%     ieegCarImpFilt(iChan,:,:) = filt60(squeeze(double(ieegCarImp(iChan,:,:))),fsD);
% end
% ieegCarImp = ieegCarImpFilt;
load('sigChannel.mat');
%%
[~,goodtrials] = remove_bad_trials(ieegCarImp,12);
goodTrialsCommon = extractCommonTrials(goodtrials(pvalsMCleanProdResp));
ieegCarClean = ieegCarImp(:,goodTrialsCommon,:);
fsDown =200;
specSamp = 50;
specDist = 20;
gInterval = 50:specDist:150;
normFactor = [];
ieegGammaBasedown = [];
for iGamma = 1:length(gInterval)
    iGamma
    [~,ieegGammaBasedown(iGamma,:,:)] = EcogExtractHighGammaTrial(double(squeeze(ieegCarBase)),fsD,fsDown,[gInterval(iGamma) gInterval(iGamma)+specSamp],[0 length(ieegCarBase)-1]./fsD,[0 length(ieegCarBase)-1]./fsD,[]);
    timeBaseDown = [0:size(ieegGammaBasedown,3)-1]/fsDown;
    normFactor(iGamma,:,:) = [mean(squeeze(ieegGammaBasedown(iGamma,:,timeBaseDown>=5&timeBaseDown<=35)),2) std(squeeze(ieegGammaBasedown(iGamma,:,timeBaseDown>=5&timeBaseDown<=35)),0,2)];

end
   timeBaseDown = [0:size(ieegGammaBasedown,3)-1]/fsDown;
  
  
ieegGamma = [];

for iTrial = 1:size(ieegCarClean,2)
    iTrial
    for iGamma = 1:length(gInterval)
      [~,ieegTmp] = EcogExtractHighGammaTrial(double(squeeze(ieegCarClean(:,iTrial,:))),fsD,fsDown,[gInterval(iGamma) gInterval(iGamma)+specSamp],tw,etwG,squeeze(normFactor(iGamma,:,:)),1);
    ieegGamma(:,iTrial,iGamma,:) = ieegTmp;    
    end
end
timeGamma = linspace(etwG(1),etwG(2),size(ieegGamma,4));
ieegGammaClean = ieegGamma;
 %%
 phonClean = phonAll(goodTrialsCommon);
   [phonUniq,~,phonId] = unique(phonClean);
phonId = phonId';

timeSelect = timeGamma>=-0.15&timeGamma<=0.15;
matSize = size(ieegGamma(pvalsMCleanProdResp,:,:,timeSelect));
ieegModel = reshape(ieegGamma(pvalsMCleanProdResp,:,:,timeSelect),[matSize(1) matSize(2) matSize(3)*matSize(4)]);
% ieegModel = cat(3,ieegGamma(:,:,timeSelect),ieegHiGamma(:,:,timeSelect),ieegHiFreq(:,:,timeSelect));
% ieegModel = ieegGammaReshape(pvalsMCleanProdResp,:,:);
            matSize = size(ieegModel);
            gammaFlat = reshape(permute(ieegModel,[2 1 3]),[matSize(2) matSize(1)*matSize(3)]);
            %gammaFlatBin = gammaFlat(goodTrialsCommon,:);
            phonClassBin = phonId;
            CmatAll = zeros(36,36);
           cvp = cvpartition(phonClassBin,'KFold',10,'Stratify',true);
            %cvp = cvpartition(phonClassBin,'LeaveOut');
            accAll = 0;
            labels = phonUniq;
            pComb = nchoosek([1:length(labels)],2);
            CmatCombAll = zeros(size(pComb,1),2,2);
            aucAll = zeros(1,length(labels));
            
            [lossVect,aucVect] = scoreSelect(gammaFlat,phonClassBin,200,1,10); % Hyper parameter tuning
                 figure;
                 plot(mean(lossVect,1))
                [~,nDim] = min(mean(lossVect,1));
                ytestAll = [];
                ypredAll = [];
            for nCv = 1:cvp.NumTestSets
                nCv
                train = cvp.training(nCv);
                test = cvp.test(nCv);
                ieegTrain = ieegModel(:,train,:);
                ieegTest = ieegModel(:,test,:);
                matTrain = size(ieegTrain);
                gTrain = reshape(permute(ieegTrain,[2 1 3]),[matTrain(2) matTrain(1)*matTrain(3)]);
                matTest = size(ieegTest);
                gTest = reshape(permute(ieegTest,[2 1 3]),[matTest(2) matTest(1)*matTest(3)]);
                %gTrain = gammaFlatBin(train,:);
                %gTest = gammaFlatBin(test,:);
                pTrain = phonClassBin(train);
                pTest = phonClassBin(test);
%                 [lossVect,aucVect] = scoreSelect(gTrain,pTrain,120,1,20); % Hyper parameter tuning
%                  figure;
%                  plot(mean(lossVect,1))
%                 [~,nDim] = min(mean(lossVect,1)); % Selecting the optimal principal components
        %        mean(squeeze(aucVect(:,nDim,:)),1)
                [lossMod,Cmat,yhat,aucVect] = pcaDecode(gTrain,gTest,pTrain,...
                               pTest,nDim);
            ytestAll = [ytestAll pTest'];
            ypredAll = [ypredAll yhat'];
            accAll = accAll + 1 - lossMod;
            %aucAll = aucAll + aucVect;
            %CmatAll = CmatAll + Cmat;
            end
            CmatAll = confusionmat(ytestAll,ypredAll);
            acc = trace(CmatAll)/sum(sum(CmatAll))
        

        figure; 
        cm = confusionchart(CmatAll,labels,'Normalization','row-normalized');
        sortClasses(cm,labels);
        set(gca,'FontSize',20);
        title(['Accuracy : ' num2str(acc*100) '%']);
        
%%

phonClean = phonAll(goodTrialsCommon);
   [phonUniq,~,phonId] = unique(phonClean);
phonId = phonId';

lolab = [ 1 2 5 6 8];
lodor = [3 4 7];
hilab = [17 20 21 28 31 32];
hidor = [18 19 29 30];
lablo = [9 10 25 33 34];
labhi = [11 12 26 27 35 36];
dorlo = [13 14 22 23];
dorhi = [15 16 24];
 phonCoArt1 = zeros(size(phonId));
 phonCoArt1(ismember(phonId,lolab)) = 1;
 phonCoArt1(ismember(phonId,lodor)) = 2;
 phonCoArt1(ismember(phonId,hilab)) = 3;
 phonCoArt1(ismember(phonId,hidor)) = 4;
 phonCoArt1(ismember(phonId,lablo)) = 5;
 phonCoArt1(ismember(phonId,labhi)) = 6;
 phonCoArt1(ismember(phonId,dorlo)) = 7;
 phonCoArt1(ismember(phonId,dorhi)) = 8;
 
 
 timeSelect = timeGamma>=-0.5&timeGamma<=0.5;
ieegGammaSelect =(ieegGamma(:,:,:,timeSelect));
ieegGammaImage = [];
for iGamma = 1:size(ieegGammaSelect,3)
    ieegGammaImage(:,:,:,:,iGamma) = getImages(squeeze(ieegGammaSelect(:,:,iGamma,:)),chanMap,1:128);
end
matSize = size(ieegGammaSelect);
 ieegModel = reshape(ieegGammaSelect,[matSize(1) matSize(2) matSize(3)*matSize(4)]);
  labels = {'low-labial','low-dorsal','high-labial','high-dorsal','labial-low','labial-high','dorsal-low','dorsal-high'};          
            phonClassBin = phonCoArt1;
            CmatAll = zeros(8,8);
for iTer = 1
%[ytestAll,ypredAll,optimDimAll] = pcaLinearRegressDecoderWrap(ieegModel,phonIndClass(:,1),[0 1],[0 1],1:200,10);
 [~,ytestAll,ypredAll,optimDimAll] = pcaLinearDecoderWrap(ieegModel,phonClassBin,[0 1],[0 1],150,10);
CmatAll = confusionmat(ytestAll,ypredAll);
acc = trace(CmatAll)/sum(sum(CmatAll))
accIter(iTer) = acc;

end


figure; 
cm = confusionchart(CmatAll,labels,'Normalization','row-normalized');
sortClasses(cm,labels);
set(gca,'FontSize',20);
title(['Accuracy : ' num2str(acc*100) '%']);

   %% Individual channel phoneme decoding
      
        accChan = [];
%sigChannel = find(pvalsMCleanProdResp);
% matSize = size(ieegGamma);
% gammaFlat = reshape(permute(ieegGamma,[2 1 3]),[matSize(2) matSize(1)*matSize(3)]);
timeSelect = timeGamma>=-0.25&timeGamma<=0.25;
matSize = size(ieegGamma(:,:,:,timeSelect));
ieegModel = reshape(ieegGamma(:,:,:,timeSelect),[matSize(1) matSize(2) matSize(3)*matSize(4)]);
CMatCat = [];
chanLab = [];
for iChan = 1:matSize(1)
    iChan
gammaFlatBin = [squeeze(ieegModel(iChan,:,:))];
phonClassBin = [phonCoArt1  ];
CmatAll = zeros(8,8);
cvp = cvpartition(phonClassBin,'KFold',15,'Stratify',true);
accAll = 0;
labels = {'low-labial','low-dorsal','high-labial','high-dorsal','labial-low','labial-high','dorsal-low','dorsal-high'};
           
% pComb = nchoosek([1:length(labels)],2);
% CmatCombAll = zeros(size(pComb,1),2,2);

[lossVect,aucVect] = scoreSelect(gammaFlatBin,phonClassBin,100,1,15); % Hyper parameter tuning
%     figure;
%     plot(mean(lossVect,1))
    [~,nDim] = min(mean(lossVect,1)); 
for nCv = 1:cvp.NumTestSets
        
    train = cvp.training(nCv);
    test = cvp.test(nCv);
    gTrain = gammaFlatBin(train,:);
    gTest = gammaFlatBin(test,:);
    pTrain = phonClassBin(train);
    pTest = phonClassBin(test);
%     [lossVect,aucVect] = scoreSelect(gTrain,pTrain,100,1,30); % Hyper parameter tuning
% %     figure;
% %     plot(mean(lossVect,1))
%     [~,nDim] = min(mean(lossVect,1)); % Selecting the optimal principal components
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

figure;
for iPhon = 1:size(accChan,2)
    subplot(4,2,iPhon)
    chanView(squeeze(accChan(:,iPhon)).*100,chanMap,selectedChannels,labels{iPhon},[12.5 100]);
    colormap(flipud(bone(4096)));
    set(gca,'FontSize',20);
    set(gca,'xtick',[])
    set(gca,'ytick',[])
   
end
%%
phonClean = phonAll(goodTrialsCommon);
   [phonUniq,~,phonId] = unique(phonClean);
phonId = phonId';

frolab = [2 5 6 17 20 21];
frodor = [3 4 18 19 ];
baclab = [1 8 28 31 32];
bacdor = [7 29 30];
labfro = [10 11 25 26 34 35];
labbac = [9 12 27 33 36];
dorfro = [14 15 23];
dorbac = [13 16 22 24];

 phonCoArt2 = zeros(size(phonId));
 phonCoArt2(ismember(phonId,frolab)) = 1;
 phonCoArt2(ismember(phonId,frodor)) = 2;
 phonCoArt2(ismember(phonId,baclab)) = 3;
 phonCoArt2(ismember(phonId,bacdor)) = 4;
 phonCoArt2(ismember(phonId,labfro)) = 5;
 phonCoArt2(ismember(phonId,labbac)) = 6;
 phonCoArt2(ismember(phonId,dorfro)) = 7;
 phonCoArt2(ismember(phonId,dorbac)) = 8;
 labels = {'front-labial','front-dorsal','back-labial','back-dorsal','labial-front','labial-back','dorsal-front','dorsal-back'};
      
 
 timeSelect = timeGamma>=-0.5&timeGamma<=0.5;
ieegGammaSelect =(ieegGamma(:,:,:,timeSelect));
ieegGammaImage = [];
for iGamma = 1:size(ieegGammaSelect,3)
    ieegGammaImage(:,:,:,:,iGamma) = getImages(squeeze(ieegGammaSelect(:,:,iGamma,:)),chanMap,1:128);
end
matSize = size(ieegGammaSelect);
 ieegModel = reshape(ieegGammaSelect,[matSize(1) matSize(2) matSize(3)*matSize(4)]);
 % labels = {'low-labial','low-dorsal','high-labial','high-dorsal','labial-low','labial-high','dorsal-low','dorsal-high'};          
            phonClassBin = phonCoArt2;
            CmatAll = zeros(8,8);
for iTer = 1
%[ytestAll,ypredAll,optimDimAll] = pcaLinearRegressDecoderWrap(ieegModel,phonIndClass(:,1),[0 1],[0 1],1:200,10);
 [~,ytestAll,ypredAll] = stmfDecodeWrap(ieegModel,phonClassBin,[0 1],[0 1],0);
CmatAll = confusionmat(ytestAll,ypredAll);
acc = trace(CmatAll)/sum(sum(CmatAll))
accIter(iTer) = acc;

end


figure; 
cm = confusionchart(CmatAll,labels,'Normalization','row-normalized');
sortClasses(cm,labels);
set(gca,'FontSize',20);
title(['Accuracy : ' num2str(acc*100) '%']);
%% Individual channel phoneme decoding
      
        accChan = [];
%sigChannel = find(pvalsMCleanProdResp);
% matSize = size(ieegGamma);
% gammaFlat = reshape(permute(ieegGamma,[2 1 3]),[matSize(2) matSize(1)*matSize(3)]);
timeSelect = timeGamma>=-0.25&timeGamma<=0.25;
matSize = size(ieegGamma(:,:,:,timeSelect));
ieegModel = reshape(ieegGamma(:,:,:,timeSelect),[matSize(1) matSize(2) matSize(3)*matSize(4)]);
CMatCat = [];
chanLab = [];
for iChan = 1:matSize(1)
    iChan
gammaFlatBin = [squeeze(ieegModel(iChan,:,:))];
phonClassBin = [phonCoArt2  ];
CmatAll = zeros(8,8);
cvp = cvpartition(phonClassBin,'KFold',15,'Stratify',true);
accAll = 0;
labels = {'front-labial','front-dorsal','back-labial','back-dorsal','labial-front','labial-back','dorsal-front','dorsal-back'};
           
% pComb = nchoosek([1:length(labels)],2);
% CmatCombAll = zeros(size(pComb,1),2,2);

[lossVect,aucVect] = scoreSelect(gammaFlatBin,phonClassBin,100,1,15); % Hyper parameter tuning
%     figure;
%     plot(mean(lossVect,1))
    [~,nDim] = min(mean(lossVect,1)); 
for nCv = 1:cvp.NumTestSets
        
    train = cvp.training(nCv);
    test = cvp.test(nCv);
    gTrain = gammaFlatBin(train,:);
    gTest = gammaFlatBin(test,:);
    pTrain = phonClassBin(train);
    pTest = phonClassBin(test);
%     [lossVect,aucVect] = scoreSelect(gTrain,pTrain,100,1,30); % Hyper parameter tuning
% %     figure;
% %     plot(mean(lossVect,1))
%     [~,nDim] = min(mean(lossVect,1)); % Selecting the optimal principal components
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

figure;
for iPhon = 1:size(accChan,2)
    subplot(4,2,iPhon)
    chanView(squeeze(accChan(:,iPhon)).*100,chanMap,selectedChannels,labels{iPhon},[12.5 100]);
    colormap(flipud(bone(4096)));
    set(gca,'FontSize',20);
    set(gca,'xtick',[])
    set(gca,'ytick',[])
   
end