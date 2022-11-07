addpath(genpath('E:\Box Sync\ECoG_Recon\matlab_code\'));
global DUKEDIR
DUKEDIR =  'E:\Box Sync\CoganLab\Data\Micro\Processed Data';
dLabels = dir(DUKEDIR);
dLabels = dLabels(3:end);

tw = [-3 3]; % time window
etwPerc = [-0.5 2.5];
etwProd = [-1.5 1.5];% epoch time window
prtw = [-0.5 0]; % preonset time window
pstw = [0.25 0.75]; % postonset time window
gammaF = [70 150]; % frequency in Hz

%% data Loading
Experiment = loadExperiment('S33');
fsD = Experiment.processing.ieeg.sample_rate;
Trials = dbTrials('S33',Experiment.recording.recording_day,'Speech_OvertMimeMove');
chanMap = Experiment.recording.channel_map;
chanMap = rot90(chanMap,2);
selectedChannels = sort(chanMap(~isnan(chanMap)))';
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
ieegAuditory = ieegAuditory(:,respId,:);
ieegSplit1 = ieegSplit1(:,respId,:);
ieegSplit2 = ieegSplit2(:,respId,:);
ieegSplit3 = ieegSplit3(:,respId,:);
phon1 = phon1(respId);
phon2 = phon2(respId);
phon3 = phon3(respId);


 ieegSplitAll = cat(2,ieegSplit1);
 phonAll = [phon1; phon2; phon3 ];
phonSequence =strcat(phon1,phon2,phon3);
trigAud = trigAuditory(respId)./fsD;
trigOns = trigOnset(respId)./fsD;
trigOns2 = trigOnset2(respId)./fsD;
trigOfs = trigOffset(respId)./fsD;

%% Common average referencing
load('impedance.mat');
higImpId = find(log10(impedance1)>6);
ieegCarImp = carFilterImpedance(ieegSplitAll,higImpId);
ieegCarBase = carFilterImpedance(ieegBase,higImpId);
ieegCarAud = carFilterImpedance(ieegAuditory,higImpId);

%% Spectrogram 
[~,goodtrialsbase] = remove_bad_trials(ieegCarAud,14);
[specCarBase]= getSpectrograms(ieegCarAud,goodtrialsbase,tw,etwProd,[1 200],prtw,pstw,[70 150],fsD,1);
 
figure;
[~,~,meanFreqChanOut] = specChanMap(specCarBase,(chanMap'),selectedChannels,[],etwProd,prtw,pstw,[1 200],[70 150],[-4 4],0,[]);
[~,goodtrials] = remove_bad_trials(ieegCarImp,14);
[specCarResp,pvalCarProd]= getSpectrograms(ieegCarImp,goodtrials,tw,etwProd,[1 200],[-1.5 -1],[-0.25 0.25],[70 150],fsD,1);
[~, pvalsCarProd] = fdr(pvalCarProd,0.05);
%     [~, pvalsMCleanProd] = fdr(pProd,0.05);
figure;
specChanMap(specCarResp,(chanMap'),selectedChannels,find(pvalsCarProd),etwProd,[-1.5 -1],[-0.25 0.25],[1 200],[70 150],[-4 4],0,meanFreqChanOut);


figure;
specChanMap(specCarResp,fliplr(chanMap),selectedChannels,find(pvalsCarProd),etwProd,[-1.5 -1],[-0.25 0.25],[1 200],[70 150],[-4 4],62,meanFreqChanOut);


%%
load('sigChannel.mat');
fsDown =200;
[~,goodtrials] = remove_bad_trials(double(ieegCarImp),14);
goodTrialsCommon = extractCommonTrials(goodtrials(pvalsMCleanProdResp));
ieegCarClean = ieegCarImp(:,goodTrialsCommon,:);
ieegCarAudClean = ieegCarAud(:,goodTrialsCommon,:);
respTime = trigOns(goodTrialsCommon)-trigAud(goodTrialsCommon);
respDuration = trigOfs(goodTrialsCommon)-trigOns(goodTrialsCommon);
ieegHiGammaBasedown = [];
ieegHiGammaBaseNorm = [];
for iTrial = 1:size(ieegCarAudClean,2)        
        [~,ieegHiGammaBasedown(:,iTrial,:)] = EcogExtractHighGammaTrial(double(squeeze(ieegCarAudClean(:,iTrial,:))),fsD,fsDown,[70 150],tw,prtw,[]);
end
for iChan = 1:size(ieegCarAudClean,1)
    normFactor5(iChan,:) = [mean2(squeeze(ieegHiGammaBasedown(iChan,:,:))) std2(squeeze(ieegHiGammaBasedown(iChan,:,:)))];
   % pChan(iChan) = permtest_sk(ieegGammaRespPower(iChan,:),ieegGammaBasePower(iChan,:),10000);
end
for iTrial = 1:size(ieegCarAudClean,2)        
        [~,ieegHiGammaBaseNorm(:,iTrial,:)] = EcogExtractHighGammaTrial(double(squeeze(ieegCarAudClean(:,iTrial,:))),fsD,fsDown,[70 150],tw,prtw,normFactor5,1);
end
%normFactor6 = [mean(ieegHiFreqBasedown(:,timeBaseDown>=5&timeBaseDown<=35),2) std(ieegHiFreqBasedown(:,timeBaseDown>=5&timeBaseDown<=35),0,2)];

ieegHiGammaProd = [];
ieegHiGammaPerc = [];
for iTrial = 1:size(ieegCarClean,2)
    iTrial
    [~,ieegTmp] = EcogExtractHighGammaTrial(double(squeeze(ieegCarAudClean(:,iTrial,:))),fsD,fsDown,[70 150],tw,etwPerc,normFactor5,1);
    ieegHiGammaPerc(:,iTrial,:) = ieegTmp;
    [~,ieegTmp] = EcogExtractHighGammaTrial(double(squeeze(ieegCarClean(:,iTrial,:))),fsD,fsDown,[70 150],tw,etwProd,normFactor5,1);
    ieegHiGammaProd(:,iTrial,:) = ieegTmp;
end
timeGammaPerc = linspace(etwPerc(1),etwPerc(2),size(ieegHiGammaPerc,3));
timeGammaProd = linspace(etwProd(1),etwProd(2),size(ieegHiGammaProd,3));
ieegHiGammaBaseNormTimeCorrect = [];

for iChan = 1:size(ieegHiGammaBaseNorm,1)
    iChan
    numTime = size(ieegHiGammaProd,3)/size(ieegHiGammaBaseNorm,3);
    for iTrial = 1:size(ieegHiGammaBaseNorm,2)
        randTrials = datasample(1:size(ieegHiGammaProd,2),numTime-1,'Replace',false);
        trials2join = squeeze(ieegHiGammaBaseNorm(iChan,randTrials,:));
        ieegHiGammaBaseNormTimeCorrect(iChan,iTrial,:) = [squeeze(ieegHiGammaBaseNorm(iChan,iTrial,:))' trials2join(:)'];
    end
end


%% Time-Series Permutation test
pChanProd =[];
pChanPerc = [];
for iChan = 1:size(ieegHiGammaProd,1)
    iChan
    [zValsRawAct, pValsRaw, actClust]=timePermCluster(squeeze(ieegHiGammaProd(iChan,:,:)),squeeze(ieegHiGammaBaseNormTimeCorrect(iChan,:,:)),1000,2,2.32);
    goodClust = find(cell2mat(actClust.Size)>actClust.perm99);
    p_masked = zeros(1,size(ieegHiGammaProd,3));
    for iClust = goodClust
        p_masked(actClust.Start{iClust}:actClust.Start{iClust}+actClust.Size{iClust}-1) = 1;
    end
    pChanProd(iChan,:) = p_masked;
    
    
    [zValsRawAct, pValsRaw, actClust]=timePermCluster(squeeze(ieegHiGammaPerc(iChan,:,:)),squeeze(ieegHiGammaBaseNormTimeCorrect(iChan,:,:)),1000,2,2.32);
    goodClust = find(cell2mat(actClust.Size)>actClust.perm99);
    p_masked = zeros(1,size(ieegHiGammaProd,3));
    for iClust = goodClust
        p_masked(actClust.Start{iClust}:actClust.Start{iClust}+actClust.Size{iClust}-1) = 1;
    end
    pChanPerc(iChan,:) = p_masked;
end
%% Perception cluster analysis
figure;
timeSeriesChanMap(pChanPerc(sigChannel,:),fliplr(chanMap),sigChannel,timeGammaPerc,[-0.5 1.5],[-0.5 2],[])

riseTimePerc = [];
widthTimePerc = [];
for iChan = 1:size(pChanPerc,1)
    riseTimeTemp = find(diff(pChanPerc(iChan,:))>0,2);
    widthTimeTemp = find(diff(pChanPerc(iChan,:))<0,2);
    if(isempty(riseTimeTemp)||isempty(widthTimeTemp))
        riseTimePerc(iChan) = nan;
        widthTimePerc(iChan) = nan;
    else
        riseTimePerc(iChan) = timeGammaPerc(riseTimeTemp(1));
        widthTimePerc(iChan) = sum(pChanPerc(iChan,timeGammaPerc>=riseTimePerc(iChan)))./fsDown;
%         if(length(widthTimeTemp)>1)
%             timeSigGap = timeGamma(riseTimeTemp(2)) - timeGamma(widthTimeTemp(1));
%             
%             if(timeSigGap<0.15)
%                 widthTime(iChan) = timeGamma(widthTimeTemp(2)) - timeGamma(riseTimeTemp(1));
%             else
%                 widthTime(iChan) = timeGamma(widthTimeTemp(1)) - timeGamma(riseTimeTemp(1));
%             end
%         else
%             widthTime(iChan) = timeGamma(widthTimeTemp(1)) - timeGamma(riseTimeTemp(1));
%         end
    end
end
riseTimePerc(riseTimePerc>1) = nan;
riseTimePerc(setdiff(selectedChannels,sigChannel)) = nan
widthTimePerc(setdiff(selectedChannels,sigChannel)) = nan


figure; 
chanView(riseTimePerc,fliplr(chanMap),selectedChannels,isnan(chanMap),[],[0 1],[],0)
cb = colorbar;
cb.Label.String = 'HG onset (s)';

set(gca,'FontSize',20);
set(gca,'xtick',[])
set(gca,'ytick',[])
axis equal
axis tight


figure; 
chanView(widthTimePerc,fliplr(chanMap),selectedChannels,isnan(chanMap),[],[0 2],[],0)
cb = colorbar;
cb.Label.String = 'HG duration (s)';

set(gca,'FontSize',20);
set(gca,'xtick',[])
set(gca,'ytick',[])
axis equal
axis tight


figure; 
scatter(riseTimePerc,widthTimePerc,'filled');
xlabel('HG onset (s)');
ylabel('HG duration (s)');
set(gca,'FontSize',15);

axis square
axis tight


%% Production cluster analysis
figure;
timeSeriesChanMap(pChanProd(sigChannel,:),fliplr(chanMap),sigChannel,timeGammaProd,[-0.5 1.5],[-1 1],[])

riseTimeProd = [];
widthTimeProd = [];
for iChan = 1:size(pChanProd,1)
    riseTimeTemp = find(diff(pChanProd(iChan,:))>0,2);
    widthTimeTemp = find(diff(pChanProd(iChan,:))<0,2);
    if(isempty(riseTimeTemp)||isempty(widthTimeTemp))
        riseTimeProd(iChan) = nan;
        widthTimeProd(iChan) = nan;
    else
        riseTimeProd(iChan) = timeGammaProd(riseTimeTemp(1));
        widthTimeProd(iChan) = sum(pChanProd(iChan,timeGammaProd>=riseTimeProd(iChan)))./fsDown;
%         if(length(widthTimeTemp)>1)
%             timeSigGap = timeGamma(riseTimeTemp(2)) - timeGamma(widthTimeTemp(1));
%             
%             if(timeSigGap<0.15)
%                 widthTime(iChan) = timeGamma(widthTimeTemp(2)) - timeGamma(riseTimeTemp(1));
%             else
%                 widthTime(iChan) = timeGamma(widthTimeTemp(1)) - timeGamma(riseTimeTemp(1));
%             end
%         else
%             widthTime(iChan) = timeGamma(widthTimeTemp(1)) - timeGamma(riseTimeTemp(1));
%         end
    end
end
riseTimeProd(setdiff(selectedChannels,sigChannel)) = nan
widthTimeProd(setdiff(selectedChannels,sigChannel)) = nan


figure; 
chanView(riseTimeProd,fliplr(chanMap),selectedChannels,isnan(chanMap),[],[-1 1],[],0)
cb = colorbar;
cb.Label.String = 'HG onset (s)';
cb.Ticks = [-1:0.25:1];
set(gca,'FontSize',20);
set(gca,'xtick',[])
set(gca,'ytick',[])
axis equal
axis tight


figure; 
chanView(widthTimeProd,fliplr(chanMap),selectedChannels,isnan(chanMap),[],[0 2],[],0)
cb = colorbar;
cb.Label.String = 'HG duration (s)';
cb.Ticks = [0:0.5:2];
set(gca,'FontSize',20);
set(gca,'xtick',[])
set(gca,'ytick',[])
axis equal
axis tight


figure; 
scatter(riseTimeProd,widthTimeProd,'filled');
xlabel('HG onset (s)');
ylabel('HG duration (s)');
set(gca,'FontSize',15);

axis square
axis tight

% save('S14_cluster_correct_sum_sig.mat','pChan','riseTime','widthTime');
%%
featTime = [riseTimePerc; widthTimePerc;  riseTimeProd; widthTimeProd]';
meanSil = [];
for iClust = 2:10
    meanSil(iClust-1) = 0;
    for iTer = 1:20
        [featId,featCentroid] = kmeans(featTime,iClust);

        silVal = silhouette(featTime,featId,'Euclidean');
        meanSil(iClust-1) =  meanSil(iClust-1) + nanmean((silVal));
    end
    meanSil(iClust-1) = meanSil(iClust-1)/20;
end

figure;
plot(2:10,meanSil,'LineWidth',1);
xlabel('K-Means Cluster');
ylabel('Silhouetter score');

[featId,featCentroid] = kmeans(featTime,2);


figure;
plot(featTime(featId==1,1),featTime(featId==1,2),'r.','MarkerSize',12)
hold on
plot(featTime(featId==2,1),featTime(featId==2,2),'b.','MarkerSize',12)

plot(featCentroid(:,1),featCentroid(:,2),'kx',...
     'MarkerSize',15,'LineWidth',3) 
legend('Cluster 1','Cluster 2','Centroids',...
       'Location','NW')
title 'Cluster Assignments and Centroids'
hold off

figure;
chanView(featId,fliplr(chanMap),selectedChannels,isnan(chanMap),'Cluster assignments',[],[],0)
colormap(fliplr(redblue));
%% Planning channel identification
maxTimePerc = [];
for iChan = 1:size(ieegHiGammaPerc,1)
    for iTrial = 1:size(ieegHiGammaPerc,2)
        maxTimePerc(iChan,iTrial) = timeGammaPerc(ieegHiGammaPerc(iChan,iTrial,:)==max(ieegHiGammaPerc(iChan,iTrial,:)));
    end
    figure;
    scatter(maxTimePerc(iChan,:),respTime);
end
micAudPath = 'E:\microECoG surgeries\S14 surgery\audioFiles';
[micClean,fsMic] = microphoneMergeTask(micAudPath,microphone);
micSplit = squeeze(splitIeeg((microphoneClean)',trigAud.*fsMic,etwPerc,fsMic));
timeMic = linspace(etwPerc(1),etwPerc(2),size(micSplit,2));
micIntensity = [];
for iTrial = 1:size(micSplit,1)
    micTrial = micSplit(iTrial,timeMic>=0&timeMic<=1);
    micIntensity(iTrial) = sqrt(sum(micTrial.*micTrial)/length(micTrial));
end



micSplitDown = resample(micSplit',fsDown,fsMic)';
micSplitGoodTrials = micSplitDown(goodTrialsCommon,:);
micIntensityGoodTrials = micIntensity(goodTrialsCommon);
timeMicDown = linspace(etwPerc(1),etwPerc(2),size(micSplitGoodTrials,2));
%% Visualizing averaged high-gamma - time series cluster correction
ieegGammaPercMean = pChanPerc;
ieegGammaProdMean = pChanProd;
etwPlot = [-1 2.5];

figure;
for iChan = 1 : size(ieegGammaPercMean,1)

        subplot(size(chanMap,1),size(chanMap,2),find(ismember(fliplr(chanMap)',selectedChannels(iChan))));
        %
        plot(timeGammaPerc(timeGammaPerc>=etwPlot(1)&timeGammaPerc<=etwPlot(2)),ieegGammaPercMean(iChan,timeGammaPerc>=etwPlot(1)&timeGammaPerc<=etwPlot(2)),'LineWidth',2);
           hold on;  
        plot(timeGammaProd(timeGammaProd>=etwPlot(1)&timeGammaProd<=etwPlot(2)),ieegGammaProdMean(iChan,timeGammaProd>=etwPlot(1)&timeGammaProd<=etwPlot(2)),'LineWidth',2);
        ylim([-0.1 1.5])
        
    
end


%% Visualizing averaged high-gamma
ieegGammaPercMean = squeeze(mean(ieegHiGammaPerc,2));
ieegGammaProdMean = squeeze(mean(ieegHiGammaProd,2));
etwPlot = [-1 1.5];

figure;
for iChan = 1 : size(ieegGammaPercMean,1)

        subplot(size(chanMap,1),size(chanMap,2),find(ismember(fliplr(chanMap)',selectedChannels(iChan))));
        %
        plot(timeGammaPerc(timeGammaPerc>=etwPlot(1)&timeGammaPerc<=etwPlot(2)),ieegGammaPercMean(iChan,timeGammaPerc>=etwPlot(1)&timeGammaPerc<=etwPlot(2)),'LineWidth',2);
           hold on;  
        plot(timeGammaProd(timeGammaProd>=etwPlot(1)&timeGammaProd<=etwPlot(2)),ieegGammaProdMean(iChan,timeGammaProd>=etwPlot(1)&timeGammaProd<=etwPlot(2)),'LineWidth',2);
        ylim([-0.1 3])
        
    
end

%% Correlation with audio signal
chanMicCorr = []
chanMicdelay = []
for iChan = 1:size(ieegHiGammaPerc,1)
    for iTrial = 1:size(ieegHiGammaPerc,2)
        [c,lags] = xcorr(micSplitGoodTrials(iTrial,:)',squeeze(ieegHiGammaPerc(iChan,iTrial,:)),0,'coeff');
        chanMicCorr(iChan,iTrial) = c;
        chanMicdelay(iChan,iTrial) = finddelay(micSplitGoodTrials(iTrial,:)',squeeze(ieegHiGammaPerc(iChan,iTrial,:)),0.25*fsDown);
    end
end

figure;
chanView(mean(chanMicCorr,2),fliplr(chanMap),selectedChannels,isnan(chanMap),'HG-Mic correlation',[0 0.5],[],0)


figure;
chanView(mean(chanMicdelay,2)./fsDown,fliplr(chanMap),selectedChannels,isnan(chanMap),'HG-Mic delay',[0 0.3],[],0)

%% NNMF factorization

%ieegGammaMean = cat(2,squeeze( mean(ieegHiGammaPerc,2)),squeeze( mean(ieegHiGammaProd,2)));
%ieegGammaMean = cat(2,pChanPerc,pChanProd);
timeSelectProd = timeGammaProd>=-1&timeGammaProd<=1;
timeSelectPerc = timeGammaPerc>=0&timeGammaPerc<=2;
timeGammaPercSelect = timeGammaPerc(timeSelectPerc);
timeGammaProdSelect = timeGammaProd(timeSelectProd);
ieegGammaMean = cat(2,squeeze(mean(ieegHiGammaPerc(:,:,timeSelectPerc),2)),squeeze(mean(ieegHiGammaProd(:,:,timeSelectProd),2)));


[n,m] = size(ieegGammaMean');

factSil = [];
for iClust = 1:10
    iClust
    factSil(iClust) = 0;
    factVar(iClust) = 0;
    for iTer = 1
        opt = statset('Maxiter',1000,'Display','final');
        [W,H,D] = nnmf(ieegGammaMean',iClust,'replicates',200,'Options',opt);
        [~,featId] = max(H);
%         silVal = silhouette(ieegGammaMean,featId,'Euclidean');
%         factSil(iClust) =  factSil(iClust) + nanmean((silVal));
        Dr=(norm(ieegGammaMean'-W*H,'fro')/sqrt(n*m)); 

        Dtot=(norm(ieegGammaMean','fro')/sqrt(n*m));
        factVar(iClust) = factVar(iClust) + 1-(Dr/Dtot).^2;
    end
%     factSil(iClust) = factSil(iClust)/1;
    factVar(iClust) = factVar(iClust)/1;
end

figure;
plot(factVar);
xlabel('NNMF factors');
ylabel('Variance Explained');


% figure;
% plot(factSil);
% xlabel('NNMF factors');
% ylabel('Silhouette score');

%%


[W,H,D] = nnmf(ieegGammaMean',2,'replicates',200);
Dr=(norm(ieegGammaMean'-W*H,'fro')/sqrt(n*m)); 
[~,featId] = max(H);
Dtot=(norm(ieegGammaMean','fro')/sqrt(n*m));
1-Dr/Dtot
figure
t = tiledlayout(1,2,'TileSpacing','compact');
bgAx = axes(t,'XTick',[],'YTick',[],'Box','off');
bgAx.Layout.TileSpan = [1 2];


ax1 = axes(t);
plot(ax1,timeGammaPercSelect,W(1:400,:)')
xline(ax1,timeGammaPercSelect(end),':');
ax1.Box = 'off';
xlim(ax1,[timeGammaPercSelect(1) timeGammaPercSelect(end)])
xlabel(ax1, 'Auditory')

ax2 = axes(t);
ax2.Layout.Tile = 2;
plot(ax2,timeGammaProdSelect,W(401:end,:)')
xline(ax2,timeGammaProdSelect(1),':');
ax2.YAxis.Visible = 'off';
ax2.Box = 'off';
xlim(ax2,[timeGammaProdSelect(1) timeGammaProdSelect(end)])
xlabel(ax2,'Production')

% Link the axes
linkaxes([ax1 ax2], 'y')





legend('Factor 1', 'Factor 2');
for iFactor = 1:size(H,1)
    figure;
chanView(H(iFactor,:)',fliplr(chanMap),selectedChannels,isnan(chanMap),['Factor : ' num2str(iFactor)],[],[],[]);
end

figure;
chanView(featId,fliplr(chanMap),selectedChannels,isnan(chanMap),'Factor assignments',[],[],0)
colormap((redblue));

figure;
chanView(timeGammaMax,fliplr(chanMap),selectedChannels,isnan(chanMap),'Max Peak time from Auditory onset',[0 1.5],[],0)


%% Phoneme Encoding
phonClean = phonAll(:,goodTrialsCommon);
phonSequenceClean = phonSequence(goodTrialsCommon);
binClass = [];
phonConsClass = [];
phonCVClass = [];
phonCvClassEnd = [];
phonIndClass = [];
phonIndClassEnd = [];
phonMannerClass = [];
phonVoiceClass = [];
vowel = "aeiou";
for iTrial = 1:size(phonClean,2)
    [binClassTemp,phonCVClassTemp,phonIndClassTemp,phonMannerClassTemp,phonVoiceClassTemp] = phonemeEncoder(phonClean{1,iTrial});
    phonIndClass(iTrial,1) = phonIndClassTemp;     
    phonCVClass(iTrial,1) = phonCVClassTemp;
    binClass(iTrial,1) = binClassTemp;
    phonVoiceClass(iTrial,1) = phonVoiceClassTemp;
    phonMannerClass(iTrial,1) = phonMannerClassTemp;
    [binClassTemp,phonCVClassTemp,phonIndClassTemp,phonMannerClassTemp,phonVoiceClassTemp] = phonemeEncoder(phonClean{2,iTrial});
    phonIndClass(iTrial,2) = phonIndClassTemp;
    phonCVClass(iTrial,2) = phonCVClassTemp;
    phonVoiceClass(iTrial,2) = phonVoiceClassTemp;
     phonMannerClass(iTrial,2) = phonMannerClassTemp;
    [binClassTemp,phonCVClassTemp,phonIndClassTemp,phonMannerClassTemp,phonVoiceClassTemp] = phonemeEncoder(phonClean{3,iTrial});
    phonIndClass(iTrial,3) = phonIndClassTemp;
    phonCVClass(iTrial,3) = phonCVClassTemp;
    phonVoiceClass(iTrial,3) = phonVoiceClassTemp;
     phonMannerClass(iTrial,3) = phonMannerClassTemp;
     
     phonid = find(strcmp([PhonemeSequencingInfoS1.TOKEN],phonSequenceClean{iTrial}));
    if(isempty(phonid))
        phontactic(iTrial,:) = nan(1,size(PhonemeSequencingInfoS1,10));
    else
        phontactic(iTrial,:) = table2array(PhonemeSequencingInfoS1(phonid,2:11));
    end
end


phonFirstSyllable = phonIndClass(:,1)';
phonFirstSyllable(phonFirstSyllable>4) = phonIndClass(phonFirstSyllable>4,2)';
%% Decoding test 

labels = {'a','ae','i','u','b','p','v','g','k'};
labels = {'b','p','v','g','k'};
labels = {'a','ae','i','u'};
 for tRange = 0.5
        timeSelect = timeGammaProd>=-1&timeGammaProd<=1;
        ieegGammaSelect =(ieegHiGammaProd(featId==1,:,timeSelect));
        
       
        %             
       
        accIter = []; phonError = [];
        CmatCat = zeros(4,4)
        for iTer = 1:5
            iTer
            [~,ytestAll,ypredAll,optimDimAll] = pcaLinearDecoderWrap(ieegGammaSelect,(phonFirstSyllable),[0 1],[0 1],[10:10:90],10,0);
           
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

        
 end
%%
phonemeUnique = unique(phonIndClass(:,1));
% timeSelect = timeGammaPerc>=0;
% chan2use = ieegHiGammaPerc(sigChannel,:,timeSelect);

chan2use = ieegHiGammaProd(sigChannel,:,:);
figure;
sgtitle('Production');
for iFeat = 1:2
    chanPhonMean = [];
for iPos = 1:3
    for iPhon = 1:length(phonemeUnique)   
        phonMean= squeeze(mean(chan2use(featId == iFeat,phonIndClass(:,iPos) == phonemeUnique(iPhon),:),2));
        chanPhonMean(iPhon,:) = phonMean(:);
    end
    Dist=pdist(chanPhonMean);
    [Y e] = cmdscale(Dist,2);


    labels = {'a','ae','i','u','b','p','v','g','k'};
    subplot(2,3,3*(iFeat-1)+iPos);
    for iPhon = 1:length(phonemeUnique)   
        scatter(Y(iPhon,1),Y(iPhon,2),'.','w');
        h=text(Y(iPhon,1),Y(iPhon,2),labels(iPhon),'FontWeight','bold');
        set(h,'FontSize',15)
        hold on;
        xlabel('MDS1');
        ylabel('MDS2');
    end
    title(['Position ' num2str(iPos)]);
    axis square
    set(gca,'FontSize',15);
    
end
end
%% MDS - time-series

phonemeUnique = unique(phonIndClass(:,1));
timeRange = 0:0.01:2;
chan2use = ieegHiGammaPerc(sigChannel,:,:);
figure;
for iFeat = 1:2
    chanPhonMean = [];
for iPos = 1:3
    ytime = [];
    for iTime = 1:length(timeRange)
        timeSelect = timeGammaPerc>=timeRange(iTime)&timeGammaPerc<=timeRange(iTime)+0.25;
        for iPhon = 1:length(phonemeUnique)   
            phonMean= squeeze(mean(chan2use(featId==iFeat,phonIndClass(:,iPos) == phonemeUnique(iPhon),timeSelect),2));
            chanPhonMean(iPhon,:) = phonMean(:);
        end
        Dist=pdist(chanPhonMean);
        [Y e] = cmdscale(Dist,2);
        ytime(iTime,:,:) = Y;
    end


    labels = {'a','ae','i','u','b','p','v','g','k'};
     subplot(2,3,3*(iFeat-1)+iPos);

    for iPhon = 1:length(phonemeUnique) 
       % figure;
        scatter(squeeze(ytime(:,iPhon,1)),squeeze(ytime(:,iPhon,2)),'.');
      %   line(squeeze(ytime(:,iPhon,1)),squeeze(ytime(:,iPhon,2)));
        h=text(squeeze(ytime(round(length(timeRange)/2),iPhon,1)),squeeze(ytime(round(length(timeRange)/2),iPhon,2)),labels(iPhon),'FontWeight','bold');
        set(h,'FontSize',15)
        hold on;
        xlabel('MDS1');
        ylabel('MDS2');
    end
    title(['Position ' num2str(iPos)]);
    axis square
    set(gca,'FontSize',20);
    
end
end

%% MDS - time-series

phonemeUnique = unique(phonIndClass(:,1));
timeRange = -1:0.01:0.75;
chan2use = ieegHiGammaProd(sigChannel,:,:);
figure;
for iFeat = 1:2
    chanPhonMean = [];
for iPos = 1:3
    ytime = [];
    for iTime = 1:length(timeRange)
        timeSelect = timeGammaProd>=timeRange(iTime)&timeGammaProd<=timeRange(iTime)+0.25;
        for iPhon = 1:length(phonemeUnique)   
            phonMean= squeeze(mean(chan2use(featId==iFeat,phonIndClass(:,iPos) == phonemeUnique(iPhon),timeSelect),2));
            chanPhonMean(iPhon,:) = phonMean(:);
        end
        Dist=pdist(chanPhonMean);
        [Y e] = cmdscale(Dist,2);
        ytime(iTime,:,:) = Y;
    end


    labels = {'a','ae','i','u','b','p','v','g','k'};
     subplot(2,3,3*(iFeat-1)+iPos);

    for iPhon = 1:length(phonemeUnique) 
       % figure;
        scatter(squeeze(ytime(:,iPhon,1)),squeeze(ytime(:,iPhon,2)),'.');
      %   line(squeeze(ytime(:,iPhon,1)),squeeze(ytime(:,iPhon,2)));
        h=text(squeeze(ytime(length(timeRange)/2,iPhon,1)),squeeze(ytime(length(timeRange)/2,iPhon,2)),labels(iPhon),'FontWeight','bold');
        set(h,'FontSize',15)
        hold on;
        xlabel('MDS1');
        ylabel('MDS2');
    end
    title(['Position ' num2str(iPos)]);
    axis square
    set(gca,'FontSize',20);
    
end
end

%% Extract audio duration

[audioDuration, audioIntensity] = extractAudioDuration(phonClean);
responseHold = respTime-audioDuration;
%% Tensor component analysis - perception production
timeSelectProd = timeGammaProd>=-1&timeGammaProd<=1;
timeSelectPerc = timeGammaPerc>=0&timeGammaPerc<=2;

ieegHiGammaProdTensor = tensor(permute(ieegHiGammaProd,[1 3 2]));
ieegHiGammaPercTensor = tensor(permute(ieegHiGammaPerc,[1 3 2]));

ieegHiGammaPercProdTensor = tensor(permute(cat(3,ieegHiGammaPerc(:,:,timeSelectPerc),ieegHiGammaProd(:,:,timeSelectProd)) ,[1 3 2]));
R = 15;
n_fits = 10;
err = [];
OPTS.tol = 1e-5;
OPTS.maxiters = 500;
factorVar= [];
similarFactor = [];
for r = 1:R
    
    for n = 1:n_fits
        % fit model
        est_factors = [];
        %ieegHiGammaPercProdTensor = tensor(permute(cat(3,ieegHiGammaPerc,ieegHiGammaProd) ,[1 3 2]));
        est_factors = cp_nmu(ieegHiGammaPercProdTensor,r,OPTS);

        % store error
        err(r,n) = norm(full(est_factors) - ieegHiGammaPercProdTensor)/norm(ieegHiGammaPercProdTensor);
       factorVar(r,n) = norm(full(est_factors))/norm(ieegHiGammaPercProdTensor);
       similarFactor(r,n) = innerprod(full(est_factors),ieegHiGammaPercProdTensor)/(norm(full(est_factors))*norm(ieegHiGammaProdTensor))
        trialFact = [];
        for iCom = 1:r
            trialFact(iCom,:) = est_factors.u{3}(:,iCom);
        end
        modelLossTemp = 0;
        for iPhon =1:3
            linearModel = crossval(fitcdiscr(trialFact',phonIndClass(:,iPhon)'));
            modelLossTemp = modelLossTemp + kfoldLoss(linearModel);
        end
        modelLoss(r,n) = modelLossTemp/3;
    end
    
end

OPTS.tol = 1e-5;
OPTS.maxiters = 500;
ieegHiGammaPercProdTensor = tensor(permute(cat(3,ieegHiGammaPerc(:,:,timeSelectPerc),ieegHiGammaProd(:,:,timeSelectProd)) ,[1 3 2]));
labels = {'a','ae','i','u','b','p','v','g','k'};

nDim =6;
est_factors = cp_nmu(ieegHiGammaPercProdTensor,nDim,OPTS);
% est_factors = cp_opt(ieegHiGammaPercProdTensor,nDim);
[elecFact,timeFact,trialFact] = viz_ktensor_update_timeSplit(est_factors,chanMap',selectedChannels,timeGammaPerc(timeSelectPerc),timeGammaProd(timeSelectProd),phonIndClass,respTime,labels,nDim);
[elecFact,timeFact,trialFact] = viz_ktensor_update_timeSplit_compile(est_factors,(chanMap'),selectedChannels,timeGammaPerc(timeSelectPerc),timeGammaProd(timeSelectProd),phonIndClass,nDim)
viz_nnmf_factors(chanMap',selectedChannels,timeGammaPercSelect,timeGammaProdSelect,featId,W',riseTimePerc,timeGammaMax)
%%
timeRange = -1:0.01:0.75;
labels = {'a','ae','i','u','b','p','v','g','k'};
accTime = zeros(3,length(timeRange))
aucTime = [];
for iTime = 1:length(timeRange)
    timeSelect = timeGammaProd>=timeRange(iTime)&timeGammaProd<=timeRange(iTime)+0.25;
    ieegGammaSelect =(ieegHiGammaProd(featId==2,:,timeSelect));

%     matSize = size(ieegGammaSelect);
%     ieegModel = reshape(ieegGammaSelect,[matSize(1) matSize(2) matSize(3)*matSize(4)]);
    %             
    for iPhone = 1:3
        iPhone
        CmatCat = zeros(9,9);
        aucCat = zeros(1,9);
        for iTer = 1:3
            [~,ytestAll,ypredAll] = pcaLinearDecoderWrap(ieegGammaSelect,phonIndClass(:,iPhone)',[0 1],[0 1],80,10,0);
        
        %[~,ytestAll,ypredAll,aucAll] = stmfDecodeWrap(ieegModel,phonIndClass(:,iPhone)',[0 1],[0 1],20,0);
            CmatAll = confusionmat(ytestAll,ypredAll);
            acc = trace(CmatAll)/sum(sum(CmatAll));        
            CmatCat = CmatCat+CmatAll;
%             aucCat = aucCat + mean(aucAll,1);
        end
        accAll = trace(CmatCat)/sum(sum(CmatCat))
        accTime(iPhone,iTime) = accAll;
%         aucTime(iPhone,iTime,:) = aucCat/iTer;
    end
    save('accTime_250ms_10ms_var80_prodChan.mat','accTime','timeRange');
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
%% Phoneme position combo
[phonemeCombo,~,phonId] = unique(phonCVClass,'rows');
ieegGammaMeanPhon = []
ieegGammaStdPhon = []
ieegGammaPhoneme = [];
for iPhoneme = 1:length(phonemeCombo)
    phonIds = find(phonId==iPhoneme);
    ieegGammaMeanTemp = squeeze(mean(ieegHiGammaProd(:,phonIds,:),2));
    %travellingWaveMovie(ieegGammaMean,fsD,chanMap,selectedChannelsClean,timeGamma,[-0.25 1],[-1.5 1.5],120,[finger{iFinger} '-clean'],'z-score');
    ieegGammaMeanPhon(iPhoneme,:) = mean(ieegGammaMeanTemp);
    ieegGammaStdPhon(iPhoneme,:) = std(ieegGammaMeanTemp,0,1);
    ieegGammaPhoneme(:,iPhoneme,:) = ieegGammaMeanTemp ;
end

tPlot = [-0.8:0.2:0.8];
phonemeComboLabel = {'low-labial-low','low-labial-high','low-dorsal-low','low-dorsal-high',...
    'high-labial-low','high-labial-high','high-dorsal-low',...
    'labial-low-labial','labial-low-dorsal','labial-high-labial','labial-high-dorsal',...
    'dorsal-low-labial','dorsal-low-dorsal','dorsal-high-labial','dorsal-high-dorsal'};
    
figure;
tiledlayout(length(phonemeCombo),length(tPlot),'TileSpacing','compact');

for iPhoneme = 1:length(phonemeCombo)
    
    for iTime = 1:length(tPlot)
        %subplot(length(phonemeCombo),length(tPlot),(iPhoneme-1)*length(tPlot)+iTime);
        nexttile;
        val2disp = squeeze(ieegGammaPhoneme(:,phonemeCombo(iPhoneme),find(timeGammaProd>=tPlot(iTime),1)));
        chanView(val2disp,(chanMap'),selectedChannels,isnan(chanMap),[],[-4 4],[],0);%[1,8,8,16]
        
         
          %ax.Position = [ax.Position(1) ax.Position(2) 0.75/length(tPlot) 0.9/length(phonemeCombo)];
        colormap((jet(4096)));
        if(iTime==1)
            ylh = ylabel((phonemeComboLabel(iPhoneme)));             
            hYLabel = get(gca,'YLabel');
            set(hYLabel,'rotation',0,'VerticalAlignment','top')
        end
        if(iPhoneme==1)
            title([num2str(round(tPlot(iTime),2)) ' s'])
        end
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        set(gca,'FontSize',10);
    end
end

%% Phoneme position combo
[phonemeCombo,~,phonId] = unique(phonCVClass,'rows');
ieegGammaMeanPhon = []
ieegGammaStdPhon = []
ieegGammaPhoneme = [];
for iPhoneme = 1:length(phonemeCombo)
    phonIds = find(phonId==iPhoneme);
    ieegGammaMeanTemp = squeeze(mean(ieegHiGammaProd(:,phonIds,:),2));
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
    
figure;
tiledlayout(4,length(tPlot),'TileSpacing','compact');

for iPhoneme = [2,7,11,12]
    
    for iTime = 1:length(tPlot)
        %subplot(length(phonemeCombo),length(tPlot),(iPhoneme-1)*length(tPlot)+iTime);
        nexttile;
        val2disp = squeeze(ieegGammaPhoneme(:,iPhoneme,find(timeGammaProd>=tPlot(iTime),1)));
        chanView(val2disp,(chanMap'),selectedChannels,isnan(chanMap),[],[-4 4],[],1);%[1,8,8,16]
        
         
          %ax.Position = [ax.Position(1) ax.Position(2) 0.75/length(tPlot) 0.9/length(phonemeCombo)];
        colormap((jet(4096)));
%         if(iTime==1)
%             ylh = ylabel((phonemeComboLabel(iPhoneme)));             
%             hYLabel = get(gca,'YLabel');
%             set(hYLabel,'rotation',0,'VerticalAlignment','top')
%         end
        if(iPhoneme==2)
            title([num2str(round(tPlot(iTime),2)) ' s'])
        end
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        set(gca,'FontSize',15);
    end
end

%% High Gamma covariance
% ieegGammaSig = abs(hilbert(ieegSplit(sigChannelPerc,goodTrialsCommon,time>=0&time<=2.75)));
elecComb = nchoosek(1:length(selectedChannels),2);
timeSelect = timeGammaProd>=-0.5&timeGammaProd<=0.5;
eucDist = [];
pitch = 1.33;
sigCov = [];
for eD = 1 : size(elecComb,1)
    eD
    [x1,y1] = find(ismember(chanMap,selectedChannels(elecComb(eD,1))));
    [x2,y2] = find(ismember(chanMap,selectedChannels(elecComb(eD,2))));
    eucDist(eD) = pitch.*sqrt((x2-x1)^2+(y2-y1)^2);
    for iPhoneme = 1:length(phonemeCombo)
        gamma1 = squeeze(ieegGammaPhoneme((elecComb(eD,1)),iPhoneme,timeSelect));
        gamma2 = squeeze(ieegGammaPhoneme((elecComb(eD,2)),iPhoneme,timeSelect));
        covtemp = cov(gamma1,gamma2);
        sigCov(eD,iPhoneme) = -2*covtemp(1,2)+(covtemp(1,1)+covtemp(2,2));
        
    end
end

 for iPhoneme = 1:length(phonemeCombo)
    eucDistUnique = unique(eucDist);
    for eDU = 1:length(eucDistUnique)
        eIDs = find(eucDist == eucDistUnique(eDU));
        sigCovMean(eDU,iPhoneme) = median(sigCov(eIDs,iPhoneme));
        sigCovStd(eDU,iPhoneme) = std(sigCov(eIDs,iPhoneme))./sqrt(length(eIDs));   
    end
 end
figure;
errorbar(eucDistUnique, sigCovMean(:,2)',sigCovStd(:,2)');
hold on;
errorbar(eucDistUnique, sigCovMean(:,7)',sigCovStd(:,7)');
errorbar(eucDistUnique, sigCovMean(:,11)',sigCovStd(:,11)');
errorbar(eucDistUnique, sigCovMean(:,12)',sigCovStd(:,12)');
axis square
xlabel('Distance (mm)');
ylabel('Semivariance');
set(gca,'FontSize',20);
xlim([0 10])
figure;
for iPhoneme = 1:length(phonemeCombo)
errorbar(eucDistUnique, sigCovMean(:,iPhoneme)',sigCovStd(:,iPhoneme)');
hold on;
end



%% High Gamma correlation
% ieegGammaSig = abs(hilbert(ieegSplit(sigChannelPerc,goodTrialsCommon,time>=0&time<=2.75)));
ieegSig = ieegHiGammaProd;
sigChannel = find(pvalsMCleanProdResp);
elecComb = nchoosek(1:length(selectedChannels),2);
eucDist = [];
pitch = 1.33;
sigCorr = [];
for eD = 1 : size(elecComb,1)
    eD
    [x1,y1] = find(ismember(chanMap,selectedChannels(elecComb(eD,1))));
    [x2,y2] = find(ismember(chanMap,selectedChannels(elecComb(eD,2))));
    eucDist(eD) = pitch.*sqrt((x2-x1)^2+(y2-y1)^2);
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
eucDistUnique = unique(eucDist);
for eDU = 1:length(eucDistUnique)
    eIDs = find(eucDist == eucDistUnique(eDU))
    sigCorrMean(eDU) = mean(sigCorr(eIDs));
    sigCorrStd(eDU) = std(sigCorr(eIDs))./sqrt(length(eIDs));    
end
save('corrHiGamma.mat','sigCorrMean','sigCorrStd','eucDistUnique');

%%

fo = fitoptions('Method','NonlinearLeastSquares');
expft = fittype('(1-a)*exp(-b*x)+a','options',fo);
[curve1,gof1] = fit(round(eucDistUnique,2)',sigCorrMean',expft)


lorenft = fittype('(1-a)*(b^2/(b^2+x^2))+a','options',fo);
[curve2,gof2] = fit(round(eucDistUnique,2)',sigCorrMean',lorenft)

xRange = 0:0.1:20;

figure;
errorbar(round(eucDistUnique,2), sigCorrMean,sigCorrStd);
% hold on;
% plot(xRange,curve1(xRange));
ylim([0 1]);
xlabel('Electrode distance (mm)');

%% Gaussian Process Factor analysis
timeSelect = timeGammaProd>=-0.5&timeGammaProd<=0.5;
% data packaging

ieegHiGammaProdDatPhon = [];
ieegHiGammaProdPhon =  abs(squeeze(ieegHiGammaProd(sigChannel,:,timeSelect)));
for iTrial = 1:size(ieegHiGammaProdPhon,2)
    ieegHiGammaProdDatPhon(iTrial).trialId = iTrial;
    ieegHiGammaProdDatPhon(iTrial).spikes = abs(squeeze(ieegHiGammaProdPhon(:,iTrial,:)));
end

runIdx = 17;
method = 'gpfa';
xDim = 2;
kernSD = 10;
binWidth = 5;


% numFolds = 5;
% for xDim = 1:10
%   neuralTraj(runIdx, ieegHiGammaProdDatPhon, 'method', method, 'xDim', xDim,... 
%                     'numFolds', numFolds);
% end
% 
% kernSD = 20; % select kernSD for two-stage methods
% plotPredErrorVsDim(runIdx, kernSD);
% 
% xDim = 5; % select state dimensionality
% plotPredErrorVsKernSD(runIdx, xDim);


result = neuralTraj(runIdx, ieegHiGammaProdDatPhon, 'method', method, 'xDim', xDim,... 
                    'kernSDList', kernSD,'binWidth',binWidth);
[estParams, seqTrain] = postprocess(result, 'kernSD', kernSD);
plotEachDimVsTime(seqTrain, 'xorth', result.binWidth);
