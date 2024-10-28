addpath(genpath('C:\Users\sd355\Box\ECoG_Recon\matlab_code\'));
global DUKEDIR
DUKEDIR = 'C:\Users\sd355\Box\CoganLab\Data\Micro\Processed Data';
dLabels = dir(DUKEDIR);
dLabels = dLabels(3:end);

tw = [-2 2]; % time window
etw = [-1.5 1.5];
etwPerc = [-0.5 2.5]
etwG = [-1.5 1.5];% epoch time window
prtw = [-0.5 0]; % preonset time window
pstw = [0.25 0.75]; % postonset time window
gammaF = [70 150]; % frequency in Hz
load('channelMap.mat');
chanMap = (chanMap');
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
trigAud = trigAuditory(respId)./fsD;
respTime  = trigOns-trigAud;
[respTimeSort,sortId] = sort(respTime);
%% Common average referencing
load('impedance.mat')
higImpId = find(log10(impedance1)>6);
ieegCarImpAll = carFilterImpedance(ieegSplitAll,higImpId);
ieegCarImp1 = carFilterImpedance(ieegSplit1,higImpId);
ieegCarImp2 = carFilterImpedance(ieegSplit2,higImpId);
ieegCarImp3 = carFilterImpedance(ieegSplit3,higImpId);
ieegCarAud = carFilterImpedance(ieegAuditory,higImpId);
clear ieegSplitAll ieegSplit1 ieegSplit2 ieegSplit3
%% Spectrogram 
[~,goodtrialsbase] = remove_bad_trials(ieegCarAud,14);
[specCarBase]= getSpectrograms(ieegCarAud,goodtrialsbase,tw,etw,[1 200],prtw,pstw,[70 150],fsD,1);
 
figure;
[~,~,meanFreqChanOut] = specChanMap(specCarBase,(chanMap),selectedChannels,[],etw,prtw,pstw,[1 200],[70 150],[0.7 1.4],0,[]);
[~,goodtrials] = remove_bad_trials(ieegCarImpAll,14);
[specCarResp,pvalCarProd]= getSpectrograms(ieegCarImpAll,goodtrials,tw,etw,[1 200],[-1.5 -1],[-0.25 0.25],[70 150],fsD,1);
[~, pvalsCarProd] = fdr(pvalCarProd,0.05);
%     [~, pvalsMCleanProd] = fdr(pProd,0.05);
figure;
[specMeanAll,specMeanPower] = specChanMap(specCarResp,(chanMap),selectedChannels,find(pvalsCarProd),etw,[-1.5 -1],[-0.25 0.25],[1 200],[70 150],[-4 4],0,meanFreqChanOut);

%% High Gamma Filtering
load('sigChannel.mat');

[~,goodtrials] = remove_bad_trials(ieegCarImpAll,14);
goodTrialsCommon = extractCommonTrials(goodtrials(pvalsMCleanProdResp));

ieegCarAudClean = ieegCarAud(:,:,:);
ieegCarImpClean1 = ieegCarImp1(:,:,:);
ieegCarImpClean2 = ieegCarImp2(:,:,:);
ieegCarImpClean3 = ieegCarImp3(:,:,:);
fsDown =200;

gInterval = [70 150];
normFactorLow = [];
normFactorGamma = [];

    ieegGammaBasedown = [];
    ieegGammaRespdown = [];
    for iTrial = 1:size(ieegCarAudClean,2)
        iTrial
       % ieegLowFrequencyBase(:,iTrial,:) = ExtractLowFrequencyWrap(double(squeeze(ieegCarAudClean(:,iTrial,:))), fsD,fsDown, [1 30] ,tw,prtw,[],[]);
        [~,ieegGammaBasedown(:,iTrial,:)] = EcogExtractHighGammaTrial(double(squeeze(ieegCarAudClean(:,iTrial,:))),fsD,fsDown,[gInterval(1) gInterval(2)],tw,prtw,[]);
        %[~,ieegGammaRespdown(:,iTrial,:)] = EcogExtractHighGammaTrial(double(squeeze(ieegCarClean(:,iTrial,:))),fsD,fsDown,[gInterval(iGamma) gInterval(iGamma)+specSamp],tw,[-0.25 0.25],[]);
    end
%     ieegGammaRespPower = squeeze(mean(ieegGammaRespdown,3));
%     ieegGammaBasePower = squeeze(mean(ieegGammaBasedown,3));
    for iChan = 1:size(ieegGammaBasedown,1)
       % normFactorLow(iChan,:) = [mean2(squeeze(ieegLowFrequencyBase(iChan,:,:))) std2(squeeze(ieegLowFrequencyBase(iChan,:,:)))];
        normFactorGamma(iChan,:) = [mean2(squeeze(ieegGammaBasedown(iChan,:,:))) std2(squeeze(ieegGammaBasedown(iChan,:,:)))];
   % pChan(iChan) = permtest_sk(ieegGammaRespPower(iChan,:),ieegGammaBasePower(iChan,:),10000);
    end
     %[p_fdr, p_masked] = fdr( pChan, 0.05);


ieegGamma = [];
ieegLowFrequency = [];
ieegHiGammaPerc = [];
for iTrial = 1:size(ieegCarImpClean1,2)
    iTrial
%     ieegLowFrequency(:,iTrial,:) = ExtractLowFrequencyWrap(double(squeeze(ieegCarImpClean1(:,iTrial,:))), fsD,fsDown, [1 30] ,tw,etwG,normFactorLow,1);
%     
    [~,ieegTmp] = EcogExtractHighGammaTrial(double(squeeze(ieegCarAudClean(:,iTrial,:))),fsD,fsDown,[70 150],tw,etwPerc,normFactorGamma,1);
    ieegHiGammaPerc(:,iTrial,:) = ieegTmp;
    [~,ieegTmp] = EcogExtractHighGammaTrial(double(squeeze(ieegCarImpClean1(:,iTrial,:))),fsD,fsDown,[gInterval(1) gInterval(2)],tw,etwG,normFactorGamma,2);
    ieegGamma(1,:,iTrial,:) = ieegTmp;   
    [~,ieegTmp] = EcogExtractHighGammaTrial(double(squeeze(ieegCarImpClean2(:,iTrial,:))),fsD,fsDown,[gInterval(1) gInterval(2)],tw,etwG,normFactorGamma,2);
    ieegGamma(2,:,iTrial,:) = ieegTmp;
    [~,ieegTmp] = EcogExtractHighGammaTrial(double(squeeze(ieegCarImpClean3(:,iTrial,:))),fsD,fsDown,[gInterval(1) gInterval(2)],tw,etwG,normFactorGamma,2);
    ieegGamma(3,:,iTrial,:) = ieegTmp;
    
end
timeGamma = linspace(etwG(1),etwG(2),size(ieegGamma,4));
timeGammaAuditory = linspace(etwPerc(1),etwPerc(2),size(ieegHiGammaPerc,3));
figure;
imagesc(timeGammaAuditory,[],squeeze(ieegHiGammaPerc(3,sortId,:)));
hold on;
scatter(respTimeSort,1:149,20,'k','filled');
caxis([0 3]);
xlabel('Time from Auditory onset')
set(gca,'YDir', 'normal');

figure;
imagesc(timeGammaAuditory,[],squeeze(ieegHiGammaPerc(112,sortId,:)));
caxis([0 3]);
hold on;
scatter(respTimeSort,1:149,20,'k','filled');
xlabel('Time from Auditory onset')
set(gca,'YDir', 'normal');

ieegTemp1 = squeeze(ieegHiGammaPerc(3,sortId,:));
ieegTemp2 = squeeze(ieegHiGammaPerc(112,sortId,:));
timeSelect = timeGammaAuditory>=0&timeGammaAuditory<=1.5;
timeGammaSelect = timeGammaAuditory(timeSelect);
for iTrial = 1:size(ieegTemp1,1)
    [~,maxId] = max(ieegTemp1(iTrial,timeSelect));
    timeMax1(iTrial) = timeGammaSelect(maxId);
    [~,maxId] = max(ieegTemp2(iTrial,timeSelect));
    timeMax2(iTrial) = timeGammaSelect(maxId);    
end
figure; histogram(timeMax1,20);
hold on;
histogram(timeMax2,20);
legend('Sensory Channel','Production Channel');
xlabel('Time from auditory onset (s)');

figure;
scatter(respTimeSort,timeMax1,'filled');
hold on
scatter(respTimeSort,timeMax2,'filled');
xlabel('Response time (s)');
ylabel('Maximum peak time');
% ieegCompile(:,:,1,:) = ieegLowFrequency;
% ieegCompile(:,:,2,:) = squeeze(ieegGamma(1,:,:,:));

%travellingWaveMovie(ieegGammaRespMean,fsD,chanMap,selectedChannels,timeGamma,[-1 0.5],[0 2],240,'S14GammaResponse','z-score');

%% Phoneme Encoding
phonClean = phonAll;
binClass = [];
phonConsClass = [];
phonCVClass = [];
phonCvClassEnd = [];
phonIndClass = [];
phonMannerClass = [];
phonIndClassEnd = [];
phontactic =[];
vowel = "aeiou";
for iTrial = 1:size(phonClean,2)
    [binClassTemp,phonCVClassTemp,phonIndClassTemp, phonMannerClassTemp] = phonemeEncoder(phonClean{1,iTrial});
    phonIndClass(iTrial,1) = phonIndClassTemp;     
    phonCVClass(iTrial,1) = phonCVClassTemp;
    phonMannerClass(iTrial,1) = phonMannerClassTemp;
    [binClassTemp,phonCVClassTemp,phonIndClassTemp] = phonemeEncoder(phonClean{2,iTrial});
    phonIndClass(iTrial,2) = phonIndClassTemp;
    phonCVClass(iTrial,2) = phonCVClassTemp;
    [binClassTemp,phonCVClassTemp,phonIndClassTemp] = phonemeEncoder(phonClean{3,iTrial});
    phonIndClass(iTrial,3) = phonIndClassTemp;
    phonCVClass(iTrial,3) = phonCVClassTemp;
%     phonid = find(strcmp([PhonemeSequencingInfoS1.TOKEN],phonSequence{iTrial}));
%     if(isempty(phonid))
%         phontactic(iTrial,:) = nan(1,size(PhonemeSequencingInfoS1,10));
%     else
%         phontactic(iTrial,:) = table2array(PhonemeSequencingInfoS1(phonid,2:11));
%     end
    
end
%% SVD - tsne - 2D
tRange = 0.25;
iPhon = 1;
nElec = [512 120 60 30];
sMeanPhoneme = [];
sRatioPhoneme = [];
sMeanClass =[];
sRatioClass =[];
for iTer = 1:50
for iSamp = 1:length(nElec)
    elecPtIds = ceil(poissonDisc2([16,8],nElec(iSamp)));
    elecPt = [];
    for iElec = 1:size(elecPtIds,1)
        elecPt(iElec) = chanMap(elecPtIds(iElec,1),elecPtIds(iElec,2));
    end
    elecPt = elecPt(~isnan(elecPt)); 
    elecPtcm = ismember(selectedChannels,elecPt);

    sum(elecPtcm)


    timeSelect = timeGamma>=-tRange&timeGamma<=tRange;
    ieegGammaSelect =squeeze(ieegGamma(1,elecPtcm,:,timeSelect));
    ieegShape = size(ieegGammaSelect);
    ieegGammaSelect = reshape(permute(ieegGammaSelect,[2 1 3]),[ieegShape(2) ieegShape(1)*ieegShape(3)]);
    [coeffIeeg,scoreIeeg,~,~,explained] = pca(ieegGammaSelect,'Centered',false);
    nModes = find(cumsum(explained)>80,1);

    Y = tsne(scoreIeeg(:,1:nModes));
    
    [s] = silhouette(Y,phonIndClass(:,iPhon)','Euclidean');
    sihouetteRatio = sum(s>0)/length(s);
    meanSilhouette = mean(s(s>0));
    sMeanPhoneme(iSamp,iTer) = meanSilhouette;
    sRatioPhoneme(iSamp,iTer) = sihouetteRatio;
    [s] = silhouette(Y,phonCVClass(:,iPhon)','Euclidean');
    sihouetteRatio = sum(s>0)/length(s);
    meanSilhouette = mean(s(s>0));
    sMeanClass(iSamp,iTer) = meanSilhouette;
    sRatioClass(iSamp,iTer) = sihouetteRatio;


%     figure(1);
%     sp = subplot(2,4,iSamp);
%     scatter(Y(:,1),Y(:,2),15,phonIndClass(:,iphon)','filled');
%     set(gca,'xtick',[])
%     set(gca,'ytick',[])
%     
%     axis square;
%     axis tight;
%     colormap(sp,turbo(9));
%     set(gca,'FontSize',15);
%     
%     
%     sp = subplot(2,4,iSamp+4);
%     scatter(Y(:,1),Y(:,2),15,phonCVClass(:,iphon)','filled');
%     set(gca,'xtick',[])
%     set(gca,'ytick',[])
%     
%     axis square;
%     axis tight;
%     colormap(sp, turbo(4));
%     set(gca,'FontSize',15);
%     
%     
%     if iSamp ==1
%         figure(2);
%         sp = subplot(1,2,1);
%         scatter(Y(:,1),Y(:,2),15,phonIndClass(:,iphon)','filled');
%         set(gca,'xtick',[])
%         set(gca,'ytick',[])
%         
%         axis square;
%         axis tight;
%         colormap(sp, turbo(9));
%         set(gca,'FontSize',15);
%         cb1 = colorbar;
%         cb1.Ticks = linspace(1.5,8.5,9);
%         labels = {'a','ae','i','u','b','p','v','g','k'};
%         cb1.TickLabels = labels;
%         set(gca,'FontSize',15);
%         sp = subplot(1,2,2);
%         scatter(Y(:,1),Y(:,2),15,phonCVClass(:,iphon)','filled');
%         set(gca,'xtick',[])
%         set(gca,'ytick',[])
%         
%         axis square;
%         axis tight;
%         colormap(sp,turbo(4));
%         set(gca,'FontSize',15);
%         cb2 = colorbar;
%         cb2.Ticks = linspace(1.5,3.5,4);
%         labels = {'low','high','labials','dorsals'};
%         cb2.TickLabels = labels;
%     end
    
end
    
end



figure;
subplot(2,1,1);
h = boxplot(sMeanPhoneme', 'symbol','');
set(h,{'linew'},{2})
axis square
set(gca,'XTickLabels', {'100','50','25','12.5'});
set(gca,'FontSize',10);
ylabel('Silhouette score');
title('Phoneme grouping');
subplot(2,1,2);
h = boxplot(sMeanClass', 'symbol','');
set(h,{'linew'},{2})
axis square
set(gca,'XTickLabels', {'100','50','25','12.5'});
set(gca,'FontSize',10);
xlabel('% Sampling');
title('Articulator grouping');
%% SVD - tsne - 3D
tRange = 0.25;

nElec = [512 120 60 30];
for iSamp = 1:length(nElec)
    elecPtIds = ceil(poissonDisc2([16,8],nElec(iSamp)));
    elecPt = [];
    for iElec = 1:size(elecPtIds,1)
        elecPt(iElec) = chanMap(elecPtIds(iElec,1),elecPtIds(iElec,2));
    end
    elecPt = elecPt(~isnan(elecPt)); 
    elecPtcm = ismember(selectedChannels,elecPt);

    sum(elecPtcm)


    timeSelect = timeGamma>=-tRange&timeGamma<=tRange;
    ieegGammaSelect =squeeze(ieegGamma(1,elecPtcm,:,timeSelect));
    ieegShape = size(ieegGammaSelect);
    ieegGammaSelect = reshape(permute(ieegGammaSelect,[2 1 3]),[ieegShape(2) ieegShape(1)*ieegShape(3)]);
    [coeffIeeg,scoreIeeg,~,~,explained] = pca(ieegGammaSelect,'Centered',false);
    nModes = find(cumsum(explained)>80,1);

    Y = tsne(scoreIeeg(:,1:nModes),'NumDimensions',3);


    figure(1);
    sp = subplot(1,4,iSamp);
    scatter3(Y(:,1),Y(:,2),Y(:,3),20,phonIndClass(:,1)','filled');
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    set(gca,'ztick',[])
    axis square;
    axis tight;
    colormap(sp,turbo(9));
    set(gca,'FontSize',15);
    
    figure(2)
    sp = subplot(1,4,iSamp);
    scatter3(Y(:,1),Y(:,2),Y(:,3),20,phonCVClass(:,1)','filled');
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    set(gca,'ztick',[])
    axis square;
    axis tight;
    colormap(sp, turbo(4));
    set(gca,'FontSize',15);
    
    
    if iSamp ==1
        figure(3);
        sp = subplot(1,2,1);
        scatter3(Y(:,1),Y(:,2),Y(:,3),20,phonIndClass(:,1)','filled');
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        set(gca,'ztick',[])
        axis square;
        axis tight;
        colormap(sp, turbo(9));
        set(gca,'FontSize',15);
        cb1 = colorbar;
        cb1.Ticks = linspace(1.5,8.5,9);
        labels = {'a','ae','i','u','b','p','v','g','k'};
        cb1.TickLabels = labels;
        set(gca,'FontSize',15);
        sp = subplot(1,2,2);
        scatter3(Y(:,1),Y(:,2),Y(:,3),20,phonCVClass(:,1)','filled');
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        set(gca,'ztick',[])
        axis square;
        axis tight;
        colormap(sp,turbo(4));
        set(gca,'FontSize',15);
        cb2 = colorbar;
        cb2.Ticks = linspace(1.5,3.5,4);
        labels = {'low','high','labials','dorsals'};
        cb2.TickLabels = labels;
    end
    
    
    
end
        %% Individual phoneme decoding
accAll = [];
iPhon = 3;
 for tRange =0.5
        timeSelect = timeGamma>=-tRange&timeGamma<=tRange;
        ieegGammaSelect =squeeze(ieegGamma(iPhon,sigChannel,:,timeSelect));
        ieegGammaAll = cat(2,squeeze(ieegGamma(1,sigChannel,:,timeSelect)),squeeze(ieegGamma(2,sigChannel,:,timeSelect)),squeeze(ieegGamma(3,sigChannel,:,timeSelect))); 
        phonIndAll = [phonIndClass(:,1)' phonIndClass(:,2)' phonIndClass(:,3)'];
        % ieegGammaImage = [];
        % for iGamma = 1:size(ieegGammaSelect,3)
        %     ieegGammaImage(:,:,:,:,iGamma) = getImages(squeeze(ieegGammaSelect(:,:,iGamma,:)),chanMap,find(pvalsMCleanProdResp));
        % end
        % ieegGammaSeries = permute(ieegGammaSelect,[2 4 1 3]);
        % save('gammaPhonemeRnn_1_seconds_v2.mat','ieegGammaSeries','ieegGammaImage','phonIndClass');
%         matSize = size(ieegGammaSelect);
%         ieegModel = reshape(ieegGammaSelect,[matSize(1) matSize(2) matSize(3)*matSize(4)]);
        %             
        labels = {'a','ae','i','u','b','p','v','g','k'};
        accIter = []; phonError = [];
        CmatCat = zeros(9,9);
        for iTer = 1:10
            iTer
            [~,ytestAll,ypredAll,optimDimAll] = pcaLinearDecoderWrap(ieegGammaSelect,phonIndClass(:,iPhon)',[0 1],[0 1],[80],20,0);
             %[~,ytestAll,ypredAll,~,stmfTemplate] = stmfDecodeWrap(ieegModel,phonIndClass(:,1)',[0 1],[0 1],20,0);
            %[~,ytestAll,ypredAll] = linearDecoder(ieegModel,phonIndClass(:,1)',[0 1],[0 1],10,0);
           
            CmatAll = confusionmat(ytestAll,ypredAll);
            CmatNorm = CmatAll./sum(CmatAll,2);
            accIter(iTer) = mean(diag(CmatNorm));
            CmatCat = CmatCat+CmatAll;
        
        end
       
        CmatCatNorm = CmatCat./sum(CmatCat,2);

        accTemp = mean(diag(CmatCatNorm));
        accUnBias = trace(CmatCat)./sum(sum(CmatCat))
        accAll = [accAll accTemp];
        figure; 
        cm = confusionchart(CmatCat,labels,'Normalization','row-normalized');
        sortClasses(cm,labels);
        set(gca,'FontSize',20);
        title(['Accuracy : ' num2str(accTemp*100) '%;' 'Window : ' num2str(tRange)]);

        [phonError,cmatvect,phonemeDistVect] = phonemeDistanceError(CmatCatNorm,1:9);
        phonError
%         cmatvect(phonemeDistVect==0) = [];
%         phonemeDistVect(phonemeDistVect==0) = [];
%         phonemeDistVect(cmatvect==0) = [];
%         cmatvect(cmatvect==0) = [];
        distModClean = fitlm(phonemeDistVect,cmatvect);


        phonemeDistVals = 1:0.5:9;

       figure
        scatter(phonemeDistVect,cmatvect,'filled')
        
        
        
        % Boot strap regression
        nErrorPoints = length(phonemeDistVect);
        nboots = 1000;
        cmatpredboot = zeros(nboots,length(phonemeDistVals));
        for iBoot = 1:nboots
            iBoot
            ix = ceil(nErrorPoints*rand(1,nErrorPoints));
            distModBoot = fitlm(phonemeDistVect(ix),cmatvect(ix));
            cmatpredboot(iBoot,:) = predict(distModBoot,phonemeDistVals');
        end
        
        modelfitP = prctile(cmatpredboot,[2.5 50 97.5],1);
        h2 = patch([phonemeDistVals fliplr(phonemeDistVals)],[modelfitP(1,:) fliplr(modelfitP(3,:))],[1 .7 .7]);
        set(h2,'EdgeColor','none');
        set(h2,'FaceAlpha',0.25);
        hold on;
        plot(phonemeDistVals',predict(distMod,phonemeDistVals'));
        xlabel('Phoneme distance (bits)')
        ylabel('Decoding error');
        %xlim([0.5 9.5]);
        set(gca,'FontSize',15);
 end
 
CmatChance = 100.*ones(9,9); 
 CmatChanceNorm = CmatChance./sum(CmatChance,2);
 [phonErrorChance] = phonemeDistanceError(CmatChanceNorm,1:9);
 %% Individual phoneme decoding: data generated
 
 

 for tRange = 0.25
        timeSelect = timeGamma>=-tRange&timeGamma<=tRange;
        ieegGammaSelect =squeeze(ieegGamma(1,:,:,timeSelect));       
        
        %             
        labels = {'a','ae','i','u','b','p','v','g','k'};
        accIter = []; phonError = [];
        CmatCat = zeros(9,9);
        cvp = cvpartition(size(phonIndClass,1), 'Holdout',0.3);
        train = cvp.training(1);
        test = cvp.test(1);
        ieegModelTrain = ieegGammaSelect(:,train,:);
        phonIndTrain = phonIndClass(train,:);
        ieegModelTest = ieegGammaSelect(:,test,:);
        phonIndTest = phonIndClass(test,:);   
        [ieegModelGen,phonIndGen] = genTraining(ieegModelTrain,phonIndTrain,200);
        
        ieegTrain = cat(2,ieegModelTrain,ieegModelGen);
        phonTrain = cat(1,phonIndTrain,phonIndGen);
        
%         ieegTrain = ieegModelTrain;
%         phonTrain = phonIndTrain;
        
        
        iPhon =2;
        for iTer = 1:10
            iTer
            [accAll,ytestAll,ypredAll,optimVarAll,aucAll] = pcaLinearDecoderWrapTrainTestSplit(ieegTrain,ieegModelTest,phonTrain(:,iPhon)',phonIndTest(:,iPhon)',[10:10:90],10,0);
            
            CmatAll = confusionmat(ytestAll,ypredAll);
            CmatNorm = CmatAll./sum(CmatAll,2)

            accIter(iTer) = mean(diag(CmatNorm));
            CmatCat = CmatCat+CmatAll;
        
        end
       
        CmatCatNorm = CmatCat./sum(CmatCat,2);

        accAll = mean(diag(CmatCatNorm));

        figure; 
        cm = confusionchart(CmatCat,labels,'Normalization','row-normalized');
        sortClasses(cm,labels);
        set(gca,'FontSize',20);
        title(['Accuracy : ' num2str(accAll*100) '%;' 'Window : ' num2str(tRange)]);

        [phonError,cmatvect,phonemeDistVect] = phonemeDistanceError(CmatCatNorm,1:9);
        phonError
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
 end
 
 
 
%% Individual phoneme decoding - time series

timeRange = -1.5:0.01:1.25;
labels = {'a','ae','i','u','b','p','v','g','k'};
accTime = zeros(3,length(timeRange))
aucTime = [];
for iTimeTrain = 1:length(timeRange)
    timeSelect = timeGamma>=timeRange(iTimeTrain)&timeGamma<=timeRange(iTimeTrain)+0.25;
    ieegGammaSelect =(ieegGamma(sigChannel(maxfactorId==2),:,1,timeSelect));

    matSize = size(ieegGammaSelect);
    ieegModel = reshape(ieegGammaSelect,[matSize(1) matSize(2) matSize(3)*matSize(4)]);
    %             
    for iPhone = 1:3
        iPhone
        CmatCat = zeros(9,9);
        aucCat = zeros(1,9);
        for iTer = 1:3
            [~,ytestAll,ypredAll] = pcaLinearDecoderWrap(ieegModel,phonIndClass(:,iPhone)',[0 1],[0 1],80,10,0);
        
        %[~,ytestAll,ypredAll,aucAll] = stmfDecodeWrap(ieegModel,phonIndClass(:,iPhone)',[0 1],[0 1],20,0);
            CmatAll = confusionmat(ytestAll,ypredAll);
            acc = trace(CmatAll)/sum(sum(CmatAll));        
            CmatCat = CmatCat+CmatAll;
%             aucCat = aucCat + mean(aucAll,1);
        end
        accAll = trace(CmatCat)/sum(sum(CmatCat))
        accTime(iPhone,iTimeTrain) = accAll;
%         aucTime(iPhone,iTime,:) = aucCat/iTer;
    end
    save('accTime_250ms_10ms_var80_planChan.mat','accTime','timeRange');
end


[accMax,accMaxId] = max(accTime');
figure;
plot(timeRange+0.125, accTime.*100,'LineWidth',2);
xlabel('Time (s)');
hold on;
xline(timeRange(accMaxId)+0.125,':','','LineWidth',2);
ylabel('Decoding Accuracy (%)');
yline(11.11, '--','chance','LineWidth',2);
set(gca,'FontSize',15);
title('Subject 1');
legend('P1','P2','P3')
ylim([0 100]);
axis square
 %% Individual phoneme decoding - Train-Test time series

timeRange = -1.5:0.025:1.25;
labels = {'a','ae','i','u','b','p','v','g','k'};
accTime = zeros(3,length(timeRange),length(timeRange));
aucTime = [];
%ieegGammaSelect = squeeze(ieegGamma(sigChannel,goodTrialsCommon,1,:));
for iTimeTrain = 1:length(timeRange)
    iTimeTrain
    tsTrain = [timeRange(iTimeTrain) timeRange(iTimeTrain)+0.25];
    for iTimeTest = 1:length(timeRange)
         tsTest = [timeRange(iTimeTest) timeRange(iTimeTest)+0.25];    

        
        %             
        for iPhone = 1:3
            iPhone
            CmatCat = zeros(9,9);
            aucCat = zeros(1,9);
            for iTer = 1
                [~,ytestAll,ypredAll] = pcaLinearDecoderWrapTrainTestSplit(squeeze(ieegGamma(iPhone,:,:,:)),squeeze(ieegGamma(1,:,:,:)),phonIndClass(:,iPhone)',etwG,tsTrain,tsTest,80,10,0);

            %[~,ytestAll,ypredAll,aucAll] = stmfDecodeWrap(ieegModel,phonIndClass(:,iPhone)',[0 1],[0 1],20,0);
                CmatAll = confusionmat(ytestAll,ypredAll);
                acc = trace(CmatAll)/sum(sum(CmatAll));        
                CmatCat = CmatCat+CmatAll;
    %             aucCat = aucCat + mean(aucAll,1);
            end
            CmatNorm = CmatCat./sum(CmatCat,2);
            accTemp = mean(diag(CmatNorm));
            accTime(iPhone,iTimeTrain,iTimeTest) = accTemp;
    %         aucTime(iPhone,iTime,:) = aucCat/iTer;
        end
        
    end
    save('accTime3D_250ms_25ms_var80_phonemevsoccurence_v3.mat','accTime','timeRange');
end
timeRangeCorrect = timeRange+0.125;
figure; subplot(1,3,1);
imagesc(timeRangeCorrect,timeRangeCorrect,squeeze(accTime(1,:,:))); 
caxis([0 0.3]);
set(gca,'YDir', 'normal');
xlim([-1 1]);
ylim([-1 1]);
axis square
set(gca,'FontSize',20);
subplot(1,3,2); imagesc(timeRangeCorrect,timeRangeCorrect,squeeze(accTime(2,:,:))); 
caxis([0 0.3]);
set(gca,'YDir', 'normal');
xlim([-1 1]);
ylim([-1 1]);
axis square
set(gca,'FontSize',20);
subplot(1,3,3); imagesc(timeRangeCorrect,timeRangeCorrect,squeeze(accTime(3,:,:))); 
caxis([0 0.3]);
set(gca,'YDir', 'normal');
xlim([-1 1]);
ylim([-1 1]);
set(gca,'FontSize',20);
axis square
%%
gaussFiltVal = 0.1;
timeRangeCorrect = timeRange+0.125;

[timeGridX,timeGridY] = meshgrid(timeRangeCorrect);
figure;

hold on;
acctemp = squeeze(accTime(1,:,:));
acctemp = imgaussfilt(acctemp,gaussFiltVal);
% maximum = max(max(acctemp));
% [y,x]=find(acctemp==maximum,1);
accMark = true(size(acctemp));

cutoff = prctile(acctemp(:),99);
accMark(acctemp<cutoff) = false;
acctempMark = acctemp;
acctempMark(acctemp<cutoff) = 0;
binaryImage = true(size(acctemp));
labeledImage = bwlabel(binaryImage);
measurements = regionprops(binaryImage, acctempMark, 'WeightedCentroid');
centerOfMass = measurements.WeightedCentroid;



% maximum = max(max(acctemp));
% [y,x]=find(acctemp==maximum,1);
[~,cont1] = contour(timeGridX,timeGridY,acctemp,[cutoff ,cutoff]);
cont1.LineWidth = 2;
cont1.LineColor = '#0072BD';
% hold on;
% plot(timeRangeCorrect(round(centerOfMass(1))), timeRangeCorrect(round(centerOfMass(2))), 'r*', 'LineWidth', 2, 'MarkerSize', 25);

xline(timeRangeCorrect(round(centerOfMass(1))),':','LineWidth',2,'Color','#0072BD');
yline(timeRangeCorrect(round(centerOfMass(2))),':','LineWidth',2,'Color','#0072BD');

acctemp = squeeze(accTime(2,:,:));
acctemp = imgaussfilt(acctemp,gaussFiltVal);
accMark = true(size(acctemp));

cutoff = prctile(acctemp(:),99);
accMark(acctemp<cutoff) = false;
% maximum = max(max(acctemp));
% [y,x]=find(acctemp==maximum,6);
acctempMark = acctemp;
acctempMark(acctemp<cutoff) = 0;

binaryImage = true(size(acctemp));
labeledImage = bwlabel(binaryImage);
measurements = regionprops(binaryImage, acctempMark, 'WeightedCentroid');
centerOfMass = measurements.WeightedCentroid;


[~,cont2] = contour(timeGridX,timeGridY,acctemp,[cutoff ,cutoff]);
cont2.LineWidth = 2;
cont2.LineColor = '#D95319';
xline(timeRangeCorrect(round(centerOfMass(1))),':','LineWidth',2,'Color','#D95319');
yline(timeRangeCorrect(round(centerOfMass(2))),':','LineWidth',2,'Color','#D95319');


acctemp = squeeze(accTime(3,:,:));
acctemp = imgaussfilt(acctemp,gaussFiltVal);
accMark = true(size(acctemp));

cutoff = prctile(acctemp(:),99);
accMark(acctemp<cutoff) = false;
acctempMark = acctemp;
acctempMark(acctemp<cutoff) = 0;
% maximum = max(max(acctemp));
% [y,x]=find(acctemp==maximum,1);
binaryImage = true(size(acctemp));
labeledImage = bwlabel(binaryImage);
measurements = regionprops(binaryImage,acctempMark, 'WeightedCentroid');
centerOfMass = measurements.WeightedCentroid;
[~,cont3] = contour(timeGridX,timeGridY,acctemp,[cutoff ,cutoff]);
cont3.LineWidth = 2;
cont3.LineColor = '#EDB120';
xline(timeRangeCorrect(round(centerOfMass(1))),':','LineWidth',2,'Color','#EDB120');
yline(timeRangeCorrect(round(centerOfMass(2))),':','LineWidth',2,'Color','#EDB120');

xlim([-1 1]);
ylim([-1 1]);
set(gca,'FontSize',20);
xlabel('Testing at utterance onset')
ylabel('Training at phoneme onset');
axis square

%% Individual phoneme decoding - Train-Test time series _ CVC and vcv

timeRange = -1.5:0.025:1.25;
labels = {'a','ae','i','u','b','p','v','g','k'};
accTimeCvc = zeros(3,length(timeRange),length(timeRange));
accTimeVcv = zeros(3,length(timeRange),length(timeRange));
aucTime = [];
cvcIds = (phonIndClass(:,1)>4);
%ieegGammaSelect = squeeze(ieegGamma(sigChannel,goodTrialsCommon,1,:));
for iTimeTrain = 1:length(timeRange)
    iTimeTrain
    tsTrain = [timeRange(iTimeTrain) timeRange(iTimeTrain)+0.25];
    for iTimeTest = 1:length(timeRange)
         tsTest = [timeRange(iTimeTest) timeRange(iTimeTest)+0.25];    
   
        for iPhone = 1:3
            iPhone
            
            for iTer = 1
                [~,ytestAll,ypredAll] = pcaLinearDecoderWrapTrainTestSplit(squeeze(ieegGamma(1,:,cvcIds,:)),squeeze(ieegGamma(iPhone,:,cvcIds,:)),phonIndClass(cvcIds,iPhone)',etwG,tsTrain,tsTest,80,10,0);

            %[~,ytestAll,ypredAll,aucAll] = stmfDecodeWrap(ieegModel,phonIndClass(:,iPhone)',[0 1],[0 1],20,0);
                CmatAll = confusionmat(ytestAll,ypredAll);
                acc = trace(CmatAll)/sum(sum(CmatAll));        
                CmatCat = CmatAll;
    %             aucCat = aucCat + mean(aucAll,1);
            end
            CmatNorm = CmatCat./sum(CmatCat,2);
            accTemp = mean(diag(CmatNorm));
            accTimeCvc(iPhone,iTimeTrain,iTimeTest) = accTemp;
            pValCvc(iPhone,iTimeTrain,iTimeTest) = StatThInv(ytestAll,accTemp*100);
    %         aucTime(iPhone,iTime,:) = aucCat/iTer;
            for iTer = 1
                [~,ytestAll,ypredAll] = pcaLinearDecoderWrapTrainTestSplit(squeeze(ieegGamma(1,:,~cvcIds,:)),squeeze(ieegGamma(iPhone,:,~cvcIds,:)),phonIndClass(~cvcIds,iPhone)',etwG,tsTrain,tsTest,80,10,0);

            %[~,ytestAll,ypredAll,aucAll] = stmfDecodeWrap(ieegModel,phonIndClass(:,iPhone)',[0 1],[0 1],20,0);
                CmatAll = confusionmat(ytestAll,ypredAll);
                acc = trace(CmatAll)/sum(sum(CmatAll));        
                CmatCat = CmatAll;
    %             aucCat = aucCat + mean(aucAll,1);
            end
            CmatNorm = CmatCat./sum(CmatCat,2);
            accTemp = mean(diag(CmatNorm));
            accTimeVcv(iPhone,iTimeTrain,iTimeTest) = accTemp;
            pValVcv(iPhone,iTimeTrain,iTimeTest) = StatThInv(ytestAll,accTemp*100);
        end
        
    end
    save('accTime3D_250ms_25ms_var80_varyingOnset_cvcv_v2.mat','accTimeCvc','accTimeVcv','pValCvc','pValVcv','timeRange');
end

timeRangeCorrect = timeRange+0.125;
figure; subplot(1,3,1);
imagesc(timeRangeCorrect,timeRangeCorrect,squeeze(pValCvc(1,:,:))); 
caxis([0 1]);
set(gca,'YDir', 'normal');
xlim([-1 1]);
ylim([-1 1]);
set(gca,'FontSize',20);
subplot(1,3,2); imagesc(timeRangeCorrect,timeRangeCorrect,squeeze(pValCvc(2,:,:))); 
caxis([0 1]);
set(gca,'YDir', 'normal');
xlim([-1 1]);
ylim([-1 1]);
set(gca,'FontSize',20);
subplot(1,3,3); imagesc(timeRangeCorrect,timeRangeCorrect,squeeze(pValCvc(3,:,:))); 
caxis([0 1]);
set(gca,'YDir', 'normal');
xlim([-1 1]);
ylim([-1 1]);
set(gca,'FontSize',20);

%%
prc4cutoff = 97.5;
gaussFiltVal = 0.5;
timeRangeCorrect = timeRange+0.125;
[timeGridX,timeGridY] = meshgrid(timeRangeCorrect);
figure;
subplot(1,2,1);
hold on;
acctemp = squeeze(accTimeCvc(1,:,:));
acctemp = imgaussfilt(acctemp,gaussFiltVal);
cutoff = prctile(acctemp(:),prc4cutoff);
maximum = max(max(acctemp));
[y,x]=find(acctemp==maximum,1);
[~,cont1] = contour(timeGridX,timeGridY,acctemp,[cutoff ,cutoff]);
cont1.LineWidth = 2;
cont1.LineColor = '#0072BD';


xline(timeRangeCorrect(x),':','LineWidth',2,'Color','#0072BD');
yline(timeRangeCorrect(y),':','LineWidth',2,'Color','#0072BD');

acctemp = squeeze(accTimeCvc(2,:,:));
acctemp = imgaussfilt(acctemp,gaussFiltVal);
cutoff = prctile(acctemp(:),prc4cutoff);
maximum = max(max(acctemp));
[y,x]=find(acctemp==maximum,1);
[~,cont2] = contour(timeGridX,timeGridY,acctemp,[cutoff ,cutoff]);
cont2.LineWidth = 2;
cont2.LineColor = '#D95319';
xline(timeRangeCorrect(x),':','LineWidth',2,'Color','#D95319');
yline(timeRangeCorrect(y),':','LineWidth',2,'Color','#D95319');


acctemp = squeeze(accTimeCvc(3,:,:));
acctemp = imgaussfilt(acctemp,gaussFiltVal);
cutoff = prctile(acctemp(:),prc4cutoff);
maximum = max(max(acctemp));
[y,x]=find(acctemp==maximum,1);
[~,cont3] = contour(timeGridX,timeGridY,acctemp,[cutoff ,cutoff]);
cont3.LineWidth = 2;
cont3.LineColor = '#EDB120';
xline(timeRangeCorrect(x),':','LineWidth',2,'Color','#EDB120');
yline(timeRangeCorrect(y),':','LineWidth',2,'Color','#EDB120');

xlim([-1 1]);
ylim([-1 1]);
set(gca,'FontSize',20);
xlabel('Testing at phoneme onset')
ylabel('Training at utterance onset');
axis square
title('CVC');

subplot(1,2,2);
hold on
timeRangeCorrect = timeRange+0.125;
[timeGridX,timeGridY] = meshgrid(timeRangeCorrect);

acctemp = squeeze(accTimeVcv(1,:,:));
acctemp = imgaussfilt(acctemp,gaussFiltVal);
cutoff = prctile(acctemp(:),prc4cutoff);
maximum = max(max(acctemp));
[y,x]=find(acctemp==maximum,1);
[~,cont1] = contour(timeGridX,timeGridY,acctemp,[cutoff ,cutoff]);
cont1.LineWidth = 2;
cont1.LineColor = '#0072BD';


xline(timeRangeCorrect(x),':','LineWidth',2,'Color','#0072BD');
yline(timeRangeCorrect(y),':','LineWidth',2,'Color','#0072BD');

acctemp = squeeze(accTimeVcv(2,:,:));
acctemp = imgaussfilt(acctemp,gaussFiltVal);
cutoff = prctile(acctemp(:),prc4cutoff);
maximum = max(max(acctemp));
[y,x]=find(acctemp==maximum,1);
[~,cont2] = contour(timeGridX,timeGridY,acctemp,[cutoff ,cutoff]);
cont2.LineWidth = 2;
cont2.LineColor = '#D95319';
xline(timeRangeCorrect(x),':','LineWidth',2,'Color','#D95319');
yline(timeRangeCorrect(y),':','LineWidth',2,'Color','#D95319');


acctemp = squeeze(accTimeVcv(3,:,:));
acctemp = imgaussfilt(acctemp,gaussFiltVal);
cutoff = prctile(acctemp(:),prc4cutoff);
maximum = max(max(acctemp));
[y,x]=find(acctemp==maximum,1);
[~,cont3] = contour(timeGridX,timeGridY,acctemp,[cutoff ,cutoff]);
cont3.LineWidth = 2;
cont3.LineColor = '#EDB120';
xline(timeRangeCorrect(x),':','LineWidth',2,'Color','#EDB120');
yline(timeRangeCorrect(y),':','LineWidth',2,'Color','#EDB120');

xlim([-1 1]);
ylim([-1 1]);
set(gca,'FontSize',20);
xlabel('Testing at phoneme onset')
ylabel('Training at utterance onset');
axis square

title('VCV');
%% Individual channel phoneme decoding - time series

timeRange = -1.5:0.01:1.25;
labels = {'a','ae','i','u','b','p','v','g','k'};
accTime = [];
aucTime = [];
for iChan = 1:size(ieegGamma,1)
    for iTimeTrain = 1:length(timeRange)
        timeSelect = timeGamma>=timeRange(iTimeTrain)&timeGamma<=timeRange(iTimeTrain)+0.25;
        ieegGammaSelect =(ieegGamma(iChan,:,:,timeSelect));

        matSize = size(ieegGammaSelect);
        ieegModel = reshape(ieegGammaSelect,[matSize(1) matSize(2) matSize(3)*matSize(4)]);
        %             
        for iPhone = 1:3

            CmatCat = zeros(4,4);
            aucCat = zeros(1,4);
            for iTer = 1:5
              [~,ytestAll,ypredAll] = pcaLinearDecoderWrap(ieegModel,phonCVClass(:,iPhone)',[0 1],[0 1],80,10,0);
               %[~,ytestAll,ypredAll] = stmfDecodeWrap(ieegModel,phonIndClass(:,iPhone)',[0 1],[0 1],20,0);
                CmatAll = confusionmat(ytestAll,ypredAll);
                acc = trace(CmatAll)/sum(sum(CmatAll));        
                CmatCat = CmatCat+CmatAll;
                %aucCat = aucCat + mean(aucAll,1);
            end
            accAll = trace(CmatCat)/sum(sum(CmatCat))
            accTime(iChan,iPhone,iTimeTrain) = accAll;
            %aucTime(iChan,iPhone,iTime,:) = aucCat/iTer;
        end
        
    end
    save('accTime_250ms_10ms_var80_indChan_articulator.mat','accTime','timeRange');
end


% figure;
% plot(timeRange+0.25, squeeze(accTime(:,3,:)),'LineWidth',2);

phonemeCombo = {'P1', 'P2', 'P3'};
tPlot = [-0.5:0.1:0.5];
figure;

for iPhoneme = 1:length(phonemeCombo)
    
    for iTimeTrain = 1:length(tPlot)
        subplot(length(phonemeCombo),length(tPlot),(iPhoneme-1)*length(tPlot)+iTimeTrain);
        
        val2disp = squeeze(accTime(:,iPhoneme,find(timeRange>=tPlot(iTimeTrain),1)));
        chanView(val2disp,(chanMap),selectedChannels,isnan(chanMap),[],[0.25 0.45],[],1);%[1,8,8,16]
         
         ax = gca;
          %ax.Position = [ax.Position(1) ax.Position(2) 0.75/length(tPlot) 0.9/length(phonemeCombo)];
        
        if(iTimeTrain==1)
            ylh = ylabel((phonemeCombo{iPhoneme}));             
            hYLabel = get(gca,'YLabel');
            set(hYLabel,'rotation',0,'VerticalAlignment','cap')
        end
        if(iPhoneme==1)
            title([num2str(round(tPlot(iTimeTrain),2)) ' s'])
        end
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        set(gca,'FontSize',15);
    end
end
%% Individual phoneme decoding - time series - vowel vs consonant

timeRange = -1.5:0.05:1;

accTimeCons = zeros(3,length(timeRange));
accTimeVow = zeros(3,length(timeRange));
for iTimeTrain = 1:length(timeRange)
    timeSelect = timeGamma>=timeRange(iTimeTrain)&timeGamma<=timeRange(iTimeTrain)+0.5;
    ieegGammaSelect =(ieegGamma(pvalsMCleanProdResp,:,:,timeSelect));

    matSize = size(ieegGammaSelect);
    ieegModel = reshape(ieegGammaSelect,[matSize(1) matSize(2) matSize(3)*matSize(4)]);
    %             
    for iPhone = 1:3
        phonIndClassTemp = phonIndClass(:,iPhone)';
        vowelid = phonIndClassTemp==1 | phonIndClassTemp==2 | phonIndClassTemp==3 | phonIndClassTemp==4;
        consid = find(phonIndClassTemp==5 | phonIndClassTemp==6 | phonIndClassTemp==7 | phonIndClassTemp==8 | phonIndClassTemp==9);

        CmatCat = zeros(4,4);
        for iTer = 1            
            %[~,ytestAll,ypredAll,optimDimAll] = pcaLinearDecoderWrap(ieegModel,phonIndClass(:,1)',[0 1],[0 1],90,5);
            [~,ytestAll,ypredAll] = stmfDecodeWrap(ieegModel(:,vowelid,:),phonIndClassTemp(vowelid),[0 1],[0 1],20,0);
            CmatAll = confusionmat(ytestAll,ypredAll);
            acc = trace(CmatAll)/sum(sum(CmatAll));        
            CmatCat = CmatCat+CmatAll;
        end
        accAll = trace(CmatCat)/sum(sum(CmatCat))
        accTimeVow(iPhone,iTimeTrain) = accAll;
        
        CmatCat = zeros(5,5);
        for iTer = 1            
            %[~,ytestAll,ypredAll,optimDimAll] = pcaLinearDecoderWrap(ieegModel,phonIndClass(:,1)',[0 1],[0 1],90,5);
            [~,ytestAll,ypredAll] = stmfDecodeWrap(ieegModel(:,consid,:),phonIndClassTemp(consid),[0 1],[0 1],20,0);
            CmatAll = confusionmat(ytestAll,ypredAll);
            acc = trace(CmatAll)/sum(sum(CmatAll));        
            CmatCat = CmatCat+CmatAll;
        end
        accAll = trace(CmatCat)/sum(sum(CmatCat))
        accTimeCons(iPhone,iTimeTrain) = accAll;
        
    end

end
%%
figure;
plot(timeRange+0.25, accTimeCons(1,:));
hold on;
plot(timeRange+0.25, accTimeVow(2,:));
hold on;
plot(timeRange+0.25, accTimeCons(3,:));
xlabel('Time (s)');
ylabel(' Accuracy');
title('Subject 1: CVC');
figure;
plot(timeRange+0.25, accTimeVow(1,:));
hold on;
plot(timeRange+0.25, accTimeCons(2,:));
hold on;
plot(timeRange+0.25, accTimeVow(3,:));
xlabel('Time (s)');
ylabel(' Accuracy');
title('Subject 1: VCV');
%% Individual phoneme decoding
        
timeSelect = timeGamma>=-0.25&timeGamma<=0.25;
ieegGammaSelect =(ieegGamma(:,:,:,timeSelect));
ieegGammaImage = [];
for iGamma = 1:size(ieegGammaSelect,3)
    ieegGammaImage(:,:,:,:,iGamma) = getImages(squeeze(ieegGammaSelect(:,:,iGamma,:)),chanMap,1:128);
end
matSize = size(ieegGammaSelect);
 ieegModel = reshape(ieegGammaSelect,[matSize(1) matSize(2) matSize(3)*matSize(4)]);
%             
labels = {'a','ae','i','u','b','p','v','g','k'};
accIter = [];
CmatCat = zeros(9,9);
for iTer = 1
%[ytestAll,ypredAll,optimDimAll] = pcaLinearRegressDecoderWrap(ieegModel,phonIndClass(:,1),[0 1],[0 1],1:200,10);
 [~,ytestAll,ypredAll] = stmfDecodeWrap(ieegModel,phonIndClass(:,1),[0 1],[0 1],0);
CmatAll = confusionmat(ytestAll,ypredAll);
acc = trace(CmatAll)/sum(sum(CmatAll))
accIter(iTer) = acc;
CmatCat = CmatCat+CmatAll;
end


figure; 
cm = confusionchart(CmatCat,labels,'Normalization','row-normalized');
sortClasses(cm,labels);
set(gca,'FontSize',20);
title(['Accuracy : ' num2str(acc*100) '%']);
         %% Individual phoneme decoding - consonant
timeSelect = timeGamma>=-0.5&timeGamma<=0.5;
vowelid = find(phonIndClass==1 | phonIndClass==2 | phonIndClass==3 | phonIndClass==4);
consid = find(phonIndClass==5 | phonIndClass==6 | phonIndClass==7 | phonIndClass==8 | phonIndClass==9);

ieegGammaSelect =(ieegGamma(:,consid,:,timeSelect));

matSize = size(ieegGammaSelect);
 ieegModel = reshape(ieegGammaSelect,[matSize(1) matSize(2) matSize(3)*matSize(4)]);
 labels = {'b','p','v','g','k'};
 phonClassBin = phonIndClass(consid);
 accIterCons = [];
CmatCatCons = zeros(5,5);
 for iTer = 1:10
%  [~,ytestAll,ypredAll] = pcaLinearDecoderWrap(ieegModel,phonClassBin,[0 1],[0 1],50,20);
[~,ytestAll,ypredAll] = stmfDecodeWrap(ieegModel,phonClassBin',[0 1],[0 1],10);
    CmatAll = confusionmat(ytestAll,ypredAll);
    acc = trace(CmatAll)/sum(sum(CmatAll))
    accIterCons(iTer) = acc
    CmatCatCons = CmatCatCons + CmatAll;
 end
%labels = {'Low','High','Labial','Dorsal'};

figure; 
cm = confusionchart(CmatCatCons,labels,'Normalization','row-normalized');
sortClasses(cm,labels);
set(gca,'FontSize',20);
title(['Accuracy : ' num2str(mean(accIterCons)*100) '%']);

% Individual phoneme decoding - vowel


ieegGammaSelect =(ieegGamma(:,vowelid,:,timeSelect));

matSize = size(ieegGammaSelect);
 ieegModel = reshape(ieegGammaSelect,[matSize(1) matSize(2) matSize(3)*matSize(4)]);
 labels = {'a','ae','i','u'};

   phonClassBin = phonIndClass(vowelid);
   accIterVow = [];
CmatCatVow = zeros(4,4);
 for iTer = 1:10
%[~,ytestAll,ypredAll] = pcaLinearDecoderWrap(ieegModel,phonClassBin,[0 1],[0 1],40,20);
[~,ytestAll,ypredAll] = stmfDecodeWrap(ieegModel,phonClassBin',[0 1],[0 1],10);
CmatAll = confusionmat(ytestAll,ypredAll);
acc = trace(CmatAll)/sum(sum(CmatAll))
accIterVow(iTer) = acc
CmatCatVow = CmatCatVow + CmatAll
%labels = {'Low','High','Labial','Dorsal'};
 end
figure; 
cm = confusionchart(CmatCatVow,labels,'Normalization','row-normalized');
sortClasses(cm,labels);
set(gca,'FontSize',20);
title(['Accuracy : ' num2str(mean(accIterVow)*100) '%']);
 %% Individual phoneme decoding - articular discriminatory
timeSelect = timeGamma>=-0.5&timeGamma<=0.5;
artimDiscId = find(phonIndClass==1 | phonIndClass==4 | phonIndClass==6| phonIndClass==9);
%consid = find(phonIndClass==5 | phonIndClass==6 | phonIndClass==7 | phonIndClass==8 | phonIndClass==9);

ieegGammaSelect =(ieegGamma(:,artimDiscId,:,timeSelect));
%ieegGammaSelect =cat(2,(ieegGamma(:,artimDiscId,:,timeSelect)), (ieegGammaBase(:,randi([1 144],1,20),:,timeSelect)));
matSize = size(ieegGammaSelect);
 ieegModel = reshape(ieegGammaSelect,[matSize(1) matSize(2) matSize(3)*matSize(4)]);
 labels = {'a','u','p','k'};
 %labels = {'a','u','b','g','rest'};
 phonClassBin = [phonIndClass(artimDiscId) ];
% phonClassBin = [phonIndClass(artimDiscId) 10.*ones(1,20)];
 accIterCons = [];
CmatCatArt = zeros(4,4);
 for iTer = 1:10
 %[~,ytestAll,ypredAll] = pcaLinearDecoderWrap(ieegModel,phonClassBin,[0 1],[0 1],50,10);
 [~,ytestAll,ypredAll] = stmfDecodeWrap(ieegModel,phonClassBin',[0 1],[0 1],10);
    CmatAll = confusionmat(ytestAll,ypredAll);
    acc = trace(CmatAll)/sum(sum(CmatAll))
    accIterCons(iTer) = acc
    CmatCatArt = CmatCatArt + CmatAll;
 end
 
 figure; 
cm = confusionchart(CmatCatArt,labels,'Normalization','row-normalized');
sortClasses(cm,labels);
set(gca,'FontSize',20);
title(['Accuracy : ' num2str(mean(accIterCons)*100) '%']);
%%
ieegSig = squeeze(ieegGamma(1,:,:,:));
elecComb = nchoosek(1:length(selectedChannels),2);
eucDist = [];
pitch = 1.33;
sigCorr = [];
for eD = 1 : size(elecComb,1)
    eD
    [x1,y1] = find(ismember(chanMap,selectedChannels(elecComb(eD,1))));
    [x2,y2] = find(ismember(chanMap,selectedChannels(elecComb(eD,2))));
    eucDist(eD) = round(pitch.*sqrt((x2-x1)^2+(y2-y1)^2),2);
    gamma1 = squeeze(ieegSig((elecComb(eD,1)),:,:));
    gamma2 = squeeze(ieegSig((elecComb(eD,2)),:,:));
    sigCorrTemp = 0;
    for iTrial = 1:size(gamma1,1)
        
        sigCorrTemp = sigCorrTemp + xcorr(gamma1(iTrial,:)',gamma2(iTrial,:)',0,'coeff');  
    end
    sigCorr(eD) = sigCorrTemp/size(gamma1,1);
end
sigCorrMean = [];
sigCorrStd = [];
sigCorrStdPos = [];
sigCorrStdNeg = [];
eucDistUnique = unique(eucDist);
for eDU = 1:length(eucDistUnique)
    eIDs = find(eucDist == eucDistUnique(eDU))
    sigCorrMean(eDU) = median(sigCorr(eIDs));
    sigCorrStd(eDU) = std(sigCorr(eIDs));   
    sigCorrStdPos(eDU) = prctile(sigCorr(eIDs),75);
     sigCorrStdNeg(eDU) = prctile(sigCorr(eIDs),25);
end



fo = fitoptions('Method','NonlinearLeastSquares');
expft = fittype('(1-a)*exp(-b*x)+a','options',fo);
[curve1,gof1] = fit(eucDist',sigCorr',expft)




figure;
scatter(eucDist,sigCorr,'filled')
hold on;
plot(eucDistUnique,curve1(eucDistUnique));
xlabel('Electrode distance (mm)')
ylabel('Correlation');
%xlim([0.5 9.5]);
set(gca,'FontSize',15);

figure;
errorbar(round(eucDistUnique,2), sigCorrMean,sigCorrMean-sigCorrStdNeg,sigCorrStdPos-sigCorrMean);
hold on;
plot(eucDistUnique,curve1(eucDistUnique));
ylim([0 1]);
xlabel('Electrode distance (mm)');

save('S14-corrHiGamma.mat','sigCorrMean','sigCorrStd','sigCorrStdNeg','sigCorrStdPos','eucDistUnique','curve1');
%% Phoneme position combo
[phonemeCombo,phonCount,phonId] = unique(phonCVClass,'rows');
ieegGammaMeanPhon = []
ieegGammaStdPhon = []
ieegGammaPhoneme = [];
for iPhoneme = 1:length(phonemeCombo)
    phonIds = find(phonId==iPhoneme);
    ieegGammaMeanTemp = squeeze(mean(ieegSig(:,phonIds,:),2));
    %travellingWaveMovie(ieegGammaMean,fsD,chanMap,selectedChannelsClean,timeGamma,[-0.25 1],[-1.5 1.5],120,[finger{iFinger} '-clean'],'z-score');
    ieegGammaMeanPhon(iPhoneme,:) = mean(ieegGammaMeanTemp);
    ieegGammaStdPhon(iPhoneme,:) = std(ieegGammaMeanTemp,0,1);
    ieegGammaPhoneme(:,iPhoneme,:) = ieegGammaMeanTemp ;
end

tPlot = [-0.6:0.2:0.6];
phonemeComboLabel = {'low-labial-low','low-labial-high','low-dorsal-low','low-dorsal-high',...
    'high-labial-low','high-labial-high','high-dorsal-low',...
    'labial-low-labial','labial-low-dorsal','labial-high-labial','labial-high-dorsal',...
    'dorsal-low-labial','dorsal-low-dorsal','dorsal-high-labial','dorsal-high-dorsal'};
phoneSelect =     find(phonCount>10);
figure;
tiledlayout(length(phoneSelect),length(tPlot),'TileSpacing','compact');

for iPhoneme = 1:length(phoneSelect)
    
    for iTime = 1:length(tPlot)
        %subplot(length(phonemeCombo),length(tPlot),(iPhoneme-1)*length(tPlot)+iTime);
        nexttile;
        val2disp = squeeze(ieegGammaPhoneme(:,phoneSelect(iPhoneme),find(timeGamma>=tPlot(iTime),1)));
        chanView(val2disp,flipud(chanMap'),selectedChannels,isnan(chanMap),[],[-4 4],[],1);%[1,8,8,16]
        
         
          %ax.Position = [ax.Position(1) ax.Position(2) 0.75/length(tPlot) 0.9/length(phonemeCombo)];
        colormap((jet(4096)));
        if(iTime==1)
            ylh = ylabel((phonemeComboLabel(phoneSelect(iPhoneme))));             
            hYLabel = get(gca,'YLabel');
            set(hYLabel,'rotation',0,'VerticalAlignment','top','HorizontalAlignment', 'right')
        end
        if(phoneSelect(iPhoneme)==1)
            title([num2str(round(tPlot(iTime),2)) ])
        end
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        set(gca,'FontSize',15);
    end
end
%% High Gamma covariance
% ieegGammaSig = abs(hilbert(ieegSplit(sigChannelPerc,goodTrialsCommon,time>=0&time<=2.75)));
elecComb = nchoosek(1:length(selectedChannels),2);
timeSelect = timeGamma>=-0.5&timeGamma<=0.5;
eucDist = [];
pitch = 1.33;
sigCov = [];
for eD = 1 : size(elecComb,1)
    eD
    [x1,y1] = find(ismember(chanMap,selectedChannels(elecComb(eD,1))));
    [x2,y2] = find(ismember(chanMap,selectedChannels(elecComb(eD,2))));
    eucDist(eD) = round(pitch.*sqrt((x2-x1)^2+(y2-y1)^2),2);
    for iPhoneme = 1:length(phonemeCombo)
        phonIds = find(phonId==iPhoneme);
        gamma1 = squeeze(ieegSig((elecComb(eD,1)),phonIds,timeSelect));
        gamma2 = squeeze(ieegSig((elecComb(eD,2)),phonIds,timeSelect));
%          covtemp = cov(gamma1,gamma2);
%          sigCov(eD,iPhoneme) = (-2*covtemp(1,2)+(covtemp(1,1)+covtemp(2,2)))/2;

        sigCorrTemp = 0;
        for iTrial = 1:size(gamma1,1)
%          covtemp = cov(gamma1(iTrial,:),gamma2(iTrial,:));
%          sigCorrTemp = sigCorrTemp + (-2*covtemp(1,2)+(covtemp(1,1)+covtemp(2,2)))/2;
            sigCorrTemp = sigCorrTemp + xcorr(gamma1(iTrial,:)',gamma2(iTrial,:)',0,'coeff');  
        end
        sigCov(eD,iPhoneme) = sigCorrTemp/size(gamma1,1);
        
    end
end

 for iPhoneme = 1:length(phonemeCombo)
    eucDistUnique = unique(eucDist);
    for eDU = 1:length(eucDistUnique)
        eIDs = find(eucDist == eucDistUnique(eDU));
        sigCovMean(eDU,iPhoneme) = median(sigCov(eIDs,iPhoneme));
        sigCovStd(eDU,iPhoneme) = std(sigCov(eIDs,iPhoneme))./sqrt(length(eIDs));   
    end
    fo = fitoptions('Method','NonlinearLeastSquares');
expft = fittype('(1-a)*exp(-b*x)+a','options',fo);
[curve1,gof1] = fit(eucDist',sigCov(:,iPhoneme),expft);
lambFit(iPhoneme) = 1/(curve1.b);
 end
 ColorSet = ColorBand(4);

figure;
errorbar(eucDistUnique, sigCovMean(:,2)',sigCovStd(:,2)','color', ColorSet(1,:));
hold on;
errorbar(eucDistUnique, sigCovMean(:,7)',sigCovStd(:,7)','color', ColorSet(2,:));
errorbar(eucDistUnique, sigCovMean(:,11)',sigCovStd(:,11)','color', ColorSet(3,:));
errorbar(eucDistUnique, sigCovMean(:,12)',sigCovStd(:,12)','color', ColorSet(4,:));
axis square
xlabel('Distance (mm)');
ylabel('Semivariance');
set(gca,'FontSize',20);

xlim([0 10])
figure;
for iPhoneme = 1:length(phoneSelect)
errorbar(eucDistUnique, sigCovMean(:,phoneSelect(iPhoneme))',sigCovStd(:,phoneSelect(iPhoneme))');
hold on;
end
ColorSet = varycolor(length(phoneSelect));
set(gca, 'ColorOrder', ColorSet);

%% 
ieegGammaMeanTemp = squeeze(mean(ieegSig,2));
tPlot = [-0.6:0.2:0.6];
figure;
tiledlayout(1,length(tPlot),'TileSpacing','compact');
ieegGammaPhoneme = [];
for iTime = 1:length(tPlot)
        %subplot(length(phonemeCombo),length(tPlot),(iPhoneme-1)*length(tPlot)+iTime);
        nexttile;
        val2disp = squeeze(ieegGammaMeanTemp(:,find(timeGamma>=tPlot(iTime),1)));
        chanView(val2disp,(chanMap),selectedChannels,isnan(chanMap),[],[-4 4],[],0.5);%[1,8,8,16]
        title([num2str(round(tPlot(iTime),2)) ])
        colormap((jet(4096)));
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        set(gca,'FontSize',15);
end
%%
for iPos = 1:3
    
    for iPhoneme = 1:9
        phonIds = find(phonIndClass(:,iPos)'==iPhoneme);
        ieegGammaMeanTemp = squeeze(mean(ieegSig(:,phonIds,:),2));
        %travellingWaveMovie(ieegGammaMean,fsD,chanMap,selectedChannelsClean,timeGamma,[-0.25 1],[-1.5 1.5],120,[finger{iFinger} '-clean'],'z-score');

        ieegGammaPhoneme(:,iPos,iPhoneme,:) = ieegGammaMeanTemp ;
    end
end

tPlot = [-0.6:0.2:0.6];
phonemeComboLabel =labels;
figure;
tiledlayout(9,length(tPlot),'TileSpacing','compact');


    
    
        for iPhoneme = 1:9
            for iTime = 1:length(tPlot)
        %subplot(length(phonemeCombo),length(tPlot),(iPhoneme-1)*length(tPlot)+iTime);
        nexttile;
        val2disp = squeeze(ieegGammaPhoneme(:,1,iPhoneme,find(timeGamma>=tPlot(iTime),1)));
        chanView(val2disp,(chanMap),selectedChannels,isnan(chanMap),[],[-4 4],[],0.5);%[1,8,8,16]
        
         
          %ax.Position = [ax.Position(1) ax.Position(2) 0.75/length(tPlot) 0.9/length(phonemeCombo)];
        colormap((jet(4096)));
        if(iTime==1)
            ylh = ylabel((phonemeComboLabel(iPhoneme)));             
            hYLabel = get(gca,'YLabel');
            set(hYLabel,'rotation',0,'VerticalAlignment','top','HorizontalAlignment', 'right')
        end
        if(iPhoneme==1)
            title([num2str(round(tPlot(iTime),2)) ])
        end
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        set(gca,'FontSize',15);
    end
    end
%%
elecComb = nchoosek(1:length(selectedChannels),2);
timeSelect = timeGamma>=-0.5&timeGamma<=0.5;
eucDist = [];
pitch = 1.33;
sigCov = [];
for eD = 1 : size(elecComb,1)
    eD
    [x1,y1] = find(ismember(chanMap,selectedChannels(elecComb(eD,1))));
    [x2,y2] = find(ismember(chanMap,selectedChannels(elecComb(eD,2))));
    eucDist(eD) = round(pitch.*sqrt((x2-x1)^2+(y2-y1)^2),2);
     for iPhoneme = 1:9
        phonIds = find(phonIndClass(:,1)'==iPhoneme);
        gamma1 = squeeze(ieegSig((elecComb(eD,1)),phonIds,timeSelect));
        gamma2 = squeeze(ieegSig((elecComb(eD,2)),phonIds,timeSelect));
%          covtemp = cov(gamma1,gamma2);
%          sigCov(eD,iPhoneme) = (-2*covtemp(1,2)+(covtemp(1,1)+covtemp(2,2)))/2;

        sigCorrTemp = 0;
        for iTrial = 1:size(gamma1,1)
%           covtemp = cov(gamma1(iTrial,:),gamma2(iTrial,:));
%           sigCorrTemp = sigCorrTemp + (-2*covtemp(1,2)+(covtemp(1,1)+covtemp(2,2)))/2;
            sigCorrTemp = sigCorrTemp + xcorr(gamma1(iTrial,:)',gamma2(iTrial,:)',0,'coeff');  
        end
        sigCov(eD,iPhoneme) = sigCorrTemp/size(gamma1,1);
        
    end
end

 for iPhoneme = 1:9
    eucDistUnique = unique(eucDist);
    for eDU = 1:length(eucDistUnique)
        eIDs = find(eucDist == eucDistUnique(eDU));
        sigCovMean(eDU,iPhoneme) = median(sigCov(eIDs,iPhoneme));
        sigCovStd(eDU,iPhoneme) = std(sigCov(eIDs,iPhoneme))./sqrt(length(eIDs));   
    end
    fo = fitoptions('Method','NonlinearLeastSquares');
%expft = fittype('(1-a)*exp(-b*x)+a','options',fo);
lorenft = fittype('(1-a)*(b^2/(b^2+x^2))+a','options',fo);
%expft = fittype('a*(1-exp(-b*x))','options',fo);
[curve1{iPhoneme},gof1] = fit(eucDist(eucDist<=10)',sigCov((eucDist<=10),iPhoneme),lorenft);
lambFit(iPhoneme) = 1/(curve1{iPhoneme}.b);
 end
  figure;
  ColorSet = distinguishable_colors(9);
% set(gca, 'ColorOrder', ColorSet);
  for iPhoneme = [1,3,5,8]
      plot(eucDistUnique,feval(curve1{iPhoneme},eucDistUnique),'Color',ColorSet(iPhoneme,:));
      hold on;
  end
 figure;
for iPhoneme = [1,3,5,8]
errorbar(eucDistUnique, sigCovMean(:,iPhoneme)',sigCovStd(:,iPhoneme)','Color',ColorSet(iPhoneme,:));
hold on;
end

%%
ieegGammaPhoneme = [];
phonIds = find(phonIndClass(:,1)'>4);
ieegGammaMeanTemp = squeeze(mean(ieegSig(:,phonIds,:),2));
%travellingWaveMovie(ieegGammaMean,fsD,chanMap,selectedChannelsClean,timeGamma,[-0.25 1],[-1.5 1.5],120,[finger{iFinger} '-clean'],'z-score');
ieegGammaPhoneme(:,1,:) = ieegGammaMeanTemp ;

phonIds = find(phonIndClass(:,1)'<=4);
ieegGammaMeanTemp = squeeze(mean(ieegSig(:,phonIds,:),2));
%travellingWaveMovie(ieegGammaMean,fsD,chanMap,selectedChannelsClean,timeGamma,[-0.25 1],[-1.5 1.5],120,[finger{iFinger} '-clean'],'z-score');
ieegGammaPhoneme(:,2,:) = ieegGammaMeanTemp ;

   

tPlot = [-0.6:0.2:0.6];
phonemeComboLabel ={'CVC','VCV'};
figure;
%tiledlayout(2,length(tPlot),'TileSpacing','compact');


    
    
        for iPhoneme = 1:2
            for iTime = 1:length(tPlot)
        subplot(length(phonemeComboLabel),length(tPlot),(iPhoneme-1)*length(tPlot)+iTime);
        %nexttile;
        val2disp = squeeze(ieegGammaPhoneme(:,iPhoneme,find(timeGamma>=tPlot(iTime),1)));
        chanView(val2disp,(chanMap),selectedChannels,isnan(chanMap),[],[-3 3],[],0.25);%[1,8,8,16]
        
         
          %ax.Position = [ax.Position(1) ax.Position(2) 0.75/length(tPlot) 0.9/length(phonemeCombo)];
        colormap((jet(4096)));
        if(iTime==1)
            ylh = ylabel((phonemeComboLabel(iPhoneme)));             
            hYLabel = get(gca,'YLabel');
            set(hYLabel,'rotation',0,'VerticalAlignment','top','HorizontalAlignment', 'right')
        end
        if(iPhoneme==1)
            title([num2str(round(tPlot(iTime),2)) ])
        end
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        set(gca,'FontSize',15);
    end
    end
%% Individual channel phoneme decoding
      
        accChan = [];
%sigChannel = find(pvalsMCleanProdResp);
% matSize = size(ieegGamma);
% gammaFlat = reshape(permute(ieegGamma,[2 1 3]),[matSize(2) matSize(1)*matSize(3)]);
timeSelect = timeGamma>=-0.5&timeGamma<=0.5;
ieegGammaSelect =squeeze(ieegGamma(1,:,:,timeSelect));
% matSize = size(ieegGamma(:,:,timeSelect));
% ieegModel = reshape(ieegGamma(:,:,timeSelect),[matSize(1) matSize(2) matSize(3)]);
CMatCat = [];
chanLabAuc = [];
chanLabAcc = [];
for iChan = 1:size(ieegGammaSelect,1)
    iChan
gammaFlatBin = [(ieegGammaSelect(iChan,:,:))];
phonClassBin = [phonIndClass(:,1)'  ];
aucAll = 0;
[~,ytestAll,ypredAll,optimDimAll,aucAll] = pcaLinearDecoderWrap(gammaFlatBin,phonClassBin,[0 1],[0 1],[80],20,0);
%[~,ytestAll,ypredAll,aucAll] = stmfDecodeWrap(gammaFlatBin,phonClassBin,[0 1],[0 1],5,1);
    
    CmatAll = confusionmat(ytestAll,ypredAll);
    CmatCatNorm = CmatAll./sum(CmatAll,2);
acc = mean(diag(CmatCatNorm));
    accChan(iChan) = acc;
    
    tp = diag(CmatCatNorm);
    chanLabAuc(iChan,:) = mean(aucAll);
    CMatCat(iChan,:,:) = CmatAll;
    chanLabAcc(iChan,:) = tp;
    
end
sigChannel = find(pvalsMCleanProdResp);
badChan = setdiff(selectedChannels,sigChannel);
badChanId = ismember(chanMap,badChan);
%labels = {'Sonorant','Voiced-plosive','Voiceless-plosive','Fricative'};
%labels = {'l','labials','dorsals'};

labels = {'low','high','labials','dorsals'};
%labels = {'a','ae','i','u','b','p','v','g','k'};
figure;
for iPhon = 1:size(chanLabAuc,2)
     subplot(2,2,iPhon)
%figure
    chanVal = chanView(squeeze(chanLabAuc(:,iPhon)),(chanMap),selectedChannels,[],labels{iPhon},[0.6 1],[],0.5);
    colormap(flipud(bone(4096)));
    chanValMark = true(size(chanVal));
    chanValMark(chanVal<0.7) = false;
    chanVal(chanVal<0.7) = 0;
    stats = regionprops(chanValMark,chanVal,'all');
   
    set(gca,'FontSize',15);
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    axis equal
    axis tight
    if(isempty(stats))
        continue;
    else
     [intensePhoneme(iPhon),maxId] = max([stats(:).MaxIntensity]+[stats(:).MeanIntensity]);
    hold on;
    scatter(stats(maxId).ConvexHull(:,1)',stats(maxId).ConvexHull(:,2)',20,'r','filled')
    hold on
    scatter(x,y,sz,'d')
    scatter(stats(maxId).Centroid(1),stats(maxId).Centroid(1),20,'x','r','filled')
    x = sqrt(stats(maxId).ConvexArea);
    phonArea = (1.33*x)^2
    end
end
% figure;
% ax1 = axes;
% [~,phonLabialObject] = chanView(squeeze(chanLabAuc(:,3)).*100,chanMap,selectedChannels,[],[],[50 100],[],0);
% redMap = interp1([0;1],[1 1 1; 1 0 0],linspace(0,1,4096));
% phonLabialObject.AlphaData = 1;
% figure;
% ax2 = axes;
% [~,phonDorsalObject] = chanView(squeeze(chanLabAuc(:,4)).*100,chanMap,selectedChannels,[],[],[50 100],[],0);
% blueMap = interp1([0;1],[1 1 1; 0 0 1],linspace(0,1,4096));
% phonDorsalObject.AlphaData = 0.5;
% linkaxes([ax1,ax2]);
% ax2.Visible = 'off'; 
% ax2.XTick = []; 
% ax2.YTick = []; 
% colormap(ax1,(redMap));
% colormap(ax2,(blueMap));
% figure;
% C = imfuse(phonLabialObject.CData,phonLabialObject.CData,'falsecolor','Scaling','joint','ColorChannels',[1 2 2]);
% imshow(C);

%% Poisson Disk Sampling - 
timeSelect = timeGamma>=-0.4&timeGamma<=0.4;
ieegGammaSelect =(ieegGamma(:,:,:,timeSelect));
%
nonSigChannel = [];
chanMapSig = chanMap;
chanMapSig(ismember(chanMap(:),nonSigChannel)) = nan;

%ieegGammaSelect = permute(ieegGammaSelect, [2 1 3 4]);
matSize = size(ieegGammaSelect);
ieegModel = reshape(ieegGammaSelect,[matSize(1) matSize(2) matSize(3)*matSize(4)]);
nSamp = 1;
firstElec = [0.025 0.05];
totalDist = 100;
accDist = [];
phonErrorDist = [];
% accDistVow = zeros(totalDist,nSamp);
numSamp = [];
for iDist = 1:totalDist
    elecSampDensity = firstElec+(iDist-1)*firstElec(1)
    iDist
    for iSamp = 1:nSamp*(totalDist+1-iDist)
        
        elecSampClose = elecSampDensity;
         nElec = round((elecSampClose(1) + (elecSampClose(2)-elecSampClose(1)) .* rand(1,1))*8*16);
            
%         
         %elecPtIds = round(poissonDisc([8,16],poisSpace,nElec));
         elecPtIds = ceil(poissonDisc2([8,16],nElec));
            

        elecPt = [];
        for iElec = 1:size(elecPtIds,1)
            elecPt(iElec) = chanMapSig(elecPtIds(iElec,1),elecPtIds(iElec,2));
        end
        elecPt = elecPt(~isnan(elecPt)); 
        elecPtcm = ismember(selectedChannels,elecPt);
        sum(elecPtcm)
%        
        CmatCat = zeros(9,9);
        for iTer = 1
            [~,ytestAll,ypredAll] = pcaLinearDecoderWrap(ieegModel(elecPt,:,:),phonIndClass(:,1)',[0 1],[0 1],80,20,0);
            %[~,ytestAll,ypredAll] = stmfDecodeWrap(ieegModel(elecPt,:,:),phonIndClassTemp,[0 1],[0 1],10,0);
            CmatAll = confusionmat(ytestAll,ypredAll);            
            CmatCat = CmatCat+CmatAll;
        end
        CmatCatNorm = CmatCat./sum(CmatCat,2);
        acc = trace(CmatCatNorm)/9;
        accDist = [accDist acc];
        phonError = phonemeDistanceError(CmatCatNorm,1:9);
        phonErrorDist = [phonErrorDist phonError];
        numSamp = [numSamp sum(elecPtcm)];
        
    end
    save('PhonemeDecodeHighResFirstPhoneme.mat','accDist', 'phonErrorDist', 'numSamp');   
    % figure; chanView(double(elecPtcm),chanMapSig,selectedChannels,isnan(chanMapSig),'Poisson sampling',[0 1],[],0);
end

% elecSampDist = [];
% elecSampStr = [];
% for iDist = 1:84
%     elecSampDist(iDist,:) = firstElec+(iDist-1)*0.025;
%     elecSampStr{iDist} = [num2str(elecSampDist(iDist,1)) '-' num2str(elecSampDist(iDist,2))];
% end


% figure; h = boxplot(accDist')
% set(h,{'linew'},{2})
% ylim([0 0.7]);
% yline(0.1111, '--','chance','LineWidth',2);
% hold on;
% set(gca,'XTick',1:2:length(elecSampStr));
% set(gca,'XTickLabel',elecSampStr(1:2:end))
% xlabel('# Electrodes');
% ylabel('decoding accuracy');
% set(gca,'FontSize',15);
% 
% figure; h = boxplot(phonErrorDist'.*17)
% set(h,{'linew'},{2})
% %ylim([0.1 0.7]);
% hold on;
% set(gca,'XTick',1:2:length(elecSampStr));
% set(gca,'XTickLabel',elecSampStr(1:2:end))
% xlabel('# Electrodes');
% ylabel('decoding error (bits)');
% set(gca,'FontSize',15);

%%
 numSampVect = round(sqrt(11*21./numSamp(:)));
 pitchUnique = unique(numSampVect);
 accMeanPitch = [];
 numPitch = [];
 accCell = [];
 for iPitch = 1:length(pitchUnique)
     accMeanPitch(iPitch) = median(accDist(numSampVect==pitchUnique(iPitch)));
     numPitch(iPitch) = sum(numSampVect==pitchUnique(iPitch));
     accCell{iPitch} = accDist(numSampVect==pitchUnique(iPitch));
 end
 
 
 
 
 h = figure('Units', 'pixels', ...
    'Position', [100 100 750 500]);
set(h, 'DefaultTextFontSize', 15);
 %scatter(numSampVect,accDist,5,'filled');
 h = boxplot(accDist.*100,numSampVect,'positions',pitchUnique', 'labels', pitchUnique', 'symbol','','Colors','k');
 set(h,{'linew'},{2})
 hold on;
 yline(max(accMeanPitch.*100).*0.95, ':','95% threshold','LineWidth',2, 'Color','r', 'LabelHorizontalAlignment','left');
 yline(0.11.*100, '--','chance','LineWidth',2, 'LabelHorizontalAlignment','left');
 plot(pitchUnique,accMeanPitch.*100,'LineWidth',2);
 ylim([0 0.7.*100]);
 xlabel('pitch (mm)');
 ylabel('P1 Accuracy (%)');
set(gca,'FontSize',20);

set(gca,'XDir','normal')
axis square

set( gca                       , ...
    'FontName'   , 'Helvetica' );


set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'      , ...
  'YGrid'       , 'off'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1         );
%% channel Map spatial sampling decoding
chanWindow = [2 4; 3 6; 4 8; 5 10; 6 12; 7 14; 8 16];
accSpaceSamp = [];
phonErrorSpaceSamp = [];
timeSelect = timeGamma>=-0.4&timeGamma<=0.4;
ieegGammaSelect =(ieegGamma(:,:,:,timeSelect));
matSize = size(ieegGammaSelect);
ieegModel = reshape(ieegGammaSelect,[matSize(1) matSize(2) matSize(3)*matSize(4)]);
for iWindow = 1:size(chanWindow,1)
    iWindow
    elecPoints = matrixSubSample(chanMapSig,chanWindow(iWindow,:));
    accSamp = [];
    phonErrorSamp = [];
    for iSamp = 1:size(elecPoints,1)
        CMatCat = zeros(9,9);
         elecPointSamp = elecPoints(iSamp,:);
        elecPointSamp = elecPointSamp(~isnan(elecPointSamp));
        for iTer = 1
            [~,ytestAll,ypredAll] = pcaLinearDecoderWrap(ieegModel(elecPointSamp,:,:),phonIndClass(:,1)',[0 1],[0 1],80,20,0);
            CmatAll = confusionmat(ytestAll,ypredAll);
            CMatCat = CMatCat + CmatAll;
        end
        CmatCatNorm = CMatCat./sum(CMatCat,2);
        acc = trace(CmatCatNorm)/9;
        
        phonError = phonemeDistanceError(CmatCatNorm,1:9);
        accSamp(iSamp) = acc; 
            
        phonErrorSamp(iSamp) = phonError;
    end
    accSpaceSamp{iWindow} = accSamp;
    phonErrorSpaceSamp{iWindow} = phonErrorSamp;
end
elecSpaceSampStr = [];
accSpaceSampAll = [];
phonErrorSpaceSampAll = [];
for iWindow = 1:size(chanWindow,1)    
    elecSpaceSampStrTemp = [num2str(chanWindow(iWindow,1)) 'x' num2str(chanWindow(iWindow,2))];
    elecSpaceSampStr = [elecSpaceSampStr; repmat({elecSpaceSampStrTemp},length(accSpaceSamp{iWindow}'),1)];
    accSpaceSampAll = [accSpaceSampAll accSpaceSamp{iWindow}];
    accSpaceMean(iWindow) = median(accSpaceSamp{iWindow});
    phonErrorSpaceSampAll = [phonErrorSpaceSampAll phonErrorSpaceSamp{iWindow}];
end 
save('PhonemeDecodeSpaceGridFirstPhoneme.mat','accSpaceSampAll','phonErrorSpaceSampAll','elecSpaceSampStr');

 h = figure('Units', 'pixels', ...
    'Position', [100 100 750 500]);
set(h, 'DefaultTextFontSize', 15); h = boxplot(accSpaceSampAll.*100,elecSpaceSampStr','symbol','','Colors','k');
set(h,{'linew'},{2})
ylim([0 0.7.*100]);
hold on;
plot(1:7,accSpaceMean.*100,'LineWidth',2);
yline(0.1111.*100, '--','chance','LineWidth',2, 'LabelHorizontalAlignment','left');
yline(max(accSpaceMean.*100).*0.95, ':','95% threshold','LineWidth',2, 'Color','r', 'LabelHorizontalAlignment','left');
xlabel('subsample grid size');
ylabel('P1 accuracy (%)');
formatTicks(gca);
%% Random trials
numTrials = 30:10:size(phonIndClass,1);
respTime = trigOfs-trigOns;
for iTrial = 1:length(numTrials)
    tRange = 0.5;
    timeSelect = timeGamma>=-tRange&timeGamma<=tRange;
    iSamp = 1;
    while(iSamp<=50)
        iSamp
        [trialSamp,trialIdx] = datasample(phonIndClass(:,1)',numTrials(iTrial),'Replace',false);
        [countTrial] = hist(trialSamp,unique(trialSamp)); %#ok<HIST>
        if(min(countTrial)<3||length(countTrial)<9)
            continue;
        end
        ieegGammaSelect =squeeze(ieegGamma(1,:,trialIdx,timeSelect));
%         matSize = size(ieegGammaSelect);
%         ieegModel = reshape(ieegGammaSelect,[matSize(1) matSize(2) matSize(3)*matSize(4)]);
        CmatCat = zeros(9,9);
        for iTer = 1
            
            [~,ytestAll,ypredAll,optimDimAll] = pcaLinearDecoderWrap(ieegGammaSelect,(phonIndClass(trialIdx,1)'),[0 1],[0 1],95,20,0);
            CmatAll = confusionmat(ytestAll,ypredAll);
            acc = trace(CmatAll)/sum(sum(CmatAll))
            CmatCat = CmatCat+CmatAll;
        end
        CmatCatNorm = CmatCat./sum(CmatCat,2);
        accTrialSamp(iTrial,iSamp) = trace(CmatCatNorm)/9;
        [phonError,cmatvect,phonemeDistVect] = phonemeDistanceError(CmatCatNorm,1:9);
        phonErrorSamp(iTrial,iSamp) = phonError;
        trainTimeSamp(iTrial,iSamp) = sum(respTime(trialIdx));
        trainSamp(iTrial,iSamp) = length(trialIdx);
        iSamp = iSamp+1;
        trainSampMean(iTrial,iSamp) = mean(countTrial);
        
    end
end


save('PhonemeDecodeNumTrialFirstPhoneme.mat','accTrialSamp','phonErrorSamp','numTrials','trainSamp');
%%
%trainTimeVect = round(trainTimeSamp(:)/60,1);
trainTimeVect = trainSamp(:);
accTrialVect = accTrialSamp(:).*100;
trainTimeUnique = unique(trainTimeVect);
 accMeanTrain = [];
 numPitch = [];
 accCell = [];
 for iTrianTime = 1:length(trainTimeUnique)
     accMeanTrain(iTrianTime) = median(accTrialVect(trainTimeVect==trainTimeUnique(iTrianTime)));
     numPitch(iTrianTime) = sum(trainTimeVect==trainTimeUnique(iTrianTime));
     accCell{iTrianTime} = accTrialVect(trainTimeVect==trainTimeUnique(iTrianTime));
 end
 

    
 h = figure;
set(h, 'DefaultTextFontSize', 10);
 boxplot(accTrialVect,trainTimeVect,'positions',trainTimeUnique', 'labels', trainTimeUnique', 'symbol','','Colors','k');
 hold on;
  plot(trainTimeUnique,accMeanTrain,'LineWidth',2);
  hold on;
 %yline((accMeanTrain(end)).*0.95, ':','95% threshold','LineWidth',2, 'Color','r');
 yline(11.11, '--','chance','LineWidth',2);
 xline(52,':','Block 1','LineWidth',2)
 xline(104,':','Block 2','LineWidth',2)
 xline(156,':','Block 3','LineWidth',2)
 xlim([25 55])
 
  ylim([0 0.7.*100]);
 xlabel('Total number of trials');
 ylabel('Decoding Accuracy (%)');
formatTicks(gca)


