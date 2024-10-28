%%

% maxIeeg = max(ieegMean');
% nonEmptyIds = maxIeeg>=0.25;
% varIeeg = var(ieegMean');
% nonEmptyIds = varIeeg>=5e-3;

ieegMeanMod = [];
ieegShuffMod = [];

% for iTrial = 1:750
%     iTrial
%     ieegMeanMod(iTrial,:,:) = smoothdata(ieegMeanSyllable(iTrial,:,:),3,'gaussian',5);
%     ieegShuffMod(iTrial,:,:) = smoothdata(ieegMeanSyllableShuffle(iTrial,:,:),3,'gaussian',5);
% end
trial2select = 1:700;
timeSelect = timeGamma>=-1;
val2include = ~isnan(ieegMeanSyllable) & ~isnan(ieegMeanSyllableShuffle);

chan2include = val2include(1,:,1);
ieegMeanMod = ieegMeanSyllable(trial2select,:,timeSelect)  - mean(ieegMeanSyllableShuffle(trial2select,:,timeSelect),1) ;
ieegShuffMod = ieegMeanSyllableShuffle(trial2select,:,timeSelect) - mean(ieegMeanSyllableShuffle(trial2select,:,timeSelect),1) ; 

% channelInfoAll = extractChannelLocation(Subject, channelNamePooled);
% chanLoc = {channelInfoAll.Location};
 %nonEmptyIds = (~cellfun('isempty', chanLoc));
%trueChan = ~contains(channelNamePooled,{'D80','D70','D103'});


% ieegMeanSyllableSelect = ieegMeanSyllable(:,:,timeSelect);
% ieegMeanSyllableShuffleSelect = ieegMeanSyllableShuffle(:,:,timeSelect);

channelNamePooledSelect = channelNamePooled;

for iTrial = 1:700
    iTrial
    ieegMeanMod(iTrial,:,:) = smoothdata(ieegMeanMod(iTrial,:,:),3,'gaussian',10);
    ieegShuffMod(iTrial,:,:) = smoothdata(ieegShuffMod(iTrial,:,:),3,'gaussian',10);
end
showStd = 1;

groupNum = 200;
timeGammaMod = timeGamma(timeSelect);
groupingIds = zeros(1,length(channelNamePooledSelect));
 %ieegMean2plot = ieegMean2plot;
timeRange = [-1 1.5];
yrange = [-0.02 0.08];
boxplace = 0.02;
boxStd = 0.005;
ylineval = -5e-3;

xlab = 'Response onset';


%%
load colors.mat
figure;
[p,n]=numSubplots(10);

hold on;
pmaskroi = [];

hemisphere = 'rh';
subplot(p(1),p(2),1)
bnLabels = extractBNlabels('ifg',hemisphere)
roiIds = extractRoiIds(Subject, channelNamePooledSelect,bnLabels);
ieegRoi = squeeze(mean(ieegMeanMod(:,roiIds,:),1));

trials2select = randperm(size(ieegMeanMod,1),groupNum);
ieegMeanRoi = squeeze(mean(mean(ieegMeanMod(:,roiIds,:),2),1))';
ieegMeanRoiTrial = squeeze(mean(ieegMeanMod(trials2select,roiIds,:),1));
[ieegMeanRoi,ieegMeanRoiSubj] = extractMeanControlSubject(ieegMeanRoiTrial,channelNamePooledSelect(roiIds));
ieegShuffleRoi = [];
for iTer = 1:1000
    trials2select = randperm(size(ieegMeanMod,1),groupNum);
   % ieegShuffleRoi(iTer, :) = squeeze(mean(mean(ieegShuffMod(trials2select,roiIds,:),2),1))';
   ieegShuffRoiTrial = squeeze(mean(ieegShuffMod(trials2select,roiIds,:),1));
    ieegShuffleRoi(iTer, :) = extractMeanControlSubject(ieegShuffRoiTrial,channelNamePooledSelect(roiIds));
end
ieegMeanSyllableRoi = squeeze(mean(ieegMeanSyllable(:,roiIds,timeSelect),2));
ieegMeanSyllableShuffleRoi = squeeze(mean(ieegMeanSyllableShuffle(:,roiIds,timeSelect),2));
%[zValsRawAct, pValsRaw, actClust] = timePermCluster(ieegMeanSyllableRoi,ieegMeanSyllableShuffleRoi,1000,1,1.645);
 [zValsRawAct, pValsRaw, actClust] = timePermClusterAfterPerm(ieegMeanRoi,ieegShuffleRoi,1,1.645);
pmask = zeros(1,size(ieegMeanRoi,2));
% [~, ~, actClust]=timePermClusterZOptim(acctimeall,acctimeshuffle,10000,1,norminv(1-pVal2Cutoff));

for iClust=1:length(actClust.Size)
    if actClust.Size{iClust}>actClust.perm95
        pmask(actClust.Start{iClust}: ...
            actClust.Start{iClust}+(actClust.Size{iClust}-1)) ...
            = 1;
    end
end
visTimePlot(timeRange,ieegMeanRoiSubj,colval = colors(4,:), labels = xlab)
hold on;
visTimePlot(timeRange,squeeze((mean(ieegShuffMod(:,roiIds,:),1))),colval = [0.5 0.5 0.5], labels = xlab)
scatter(timeGammaMod(logical(pmask)),ylineval.*ones(1,sum(pmask)),'filled',MarkerEdgeColor=colors(4,:),MarkerFaceColor=colors(4,:));

groupingIds(roiIds) = 4;
title(['IFG ' num2str(sum(roiIds))])
ieegRoiEpochScaled = ieegRoi.*pmask;
[maxAcc,maxId] = max(ieegRoiEpochScaled');
time_centroid = [];
for iElec = 1:size(ieegRoiEpochScaled,1)
    time_centroid = [time_centroid calculateTimeCentroid(timeGammaMod,ieegRoiEpochScaled(iElec,:))];
   
end
time_centroid = time_centroid(time_centroid>-0.9);
timeSyllableCentroid{1} = time_centroid;
timeSyllableMax{1} = timeGammaMod(maxId);

if(showBoxPlot&&sum(pmask)>0)
    boxchart(boxplace*ones(size(time_centroid)),...
        double(time_centroid),'orientation','horizontal','BoxFaceColor',...
        colors(4,:)*0.75,'BoxWidth',boxStd, 'BoxLineColor',...
        colors(4,:),'MarkerStyle','none')
end
ylim(yrange)
pmaskroi(1,:) = pmask;


subplot(p(1),p(2),2)
bnLabels = extractBNlabels('rmfg',hemisphere)
roiIds = extractRoiIds(Subject, channelNamePooledSelect,bnLabels);
ieegRoi = squeeze(mean(ieegMeanMod(:,roiIds,:),1));

trials2select = randperm(size(ieegMeanMod,1),groupNum);
ieegMeanRoi = squeeze(mean(mean(ieegMeanMod(:,roiIds,:),2),1))';
ieegMeanRoiTrial = squeeze(mean(ieegMeanMod(trials2select,roiIds,:),1));
[ieegMeanRoi,ieegMeanRoiSubj] = extractMeanControlSubject(ieegMeanRoiTrial,channelNamePooledSelect(roiIds));
ieegShuffleRoi = [];
for iTer = 1:1000
    trials2select = randperm(size(ieegMeanMod,1),groupNum);
    % ieegShuffleRoi(iTer, :) = squeeze(mean(mean(ieegShuffMod(trials2select,roiIds,:),2),1))';
   ieegShuffRoiTrial = squeeze(mean(ieegShuffMod(trials2select,roiIds,:),1));
    ieegShuffleRoi(iTer, :) = extractMeanControlSubject(ieegShuffRoiTrial,channelNamePooledSelect(roiIds));
end
ieegMeanSyllableRoi = squeeze(mean(ieegMeanSyllable(:,roiIds,timeSelect),2));
ieegMeanSyllableShuffleRoi = squeeze(mean(ieegMeanSyllableShuffle(:,roiIds,timeSelect),2));
% [zValsRawAct, pValsRaw, actClust] = timePermCluster(ieegMeanSyllableRoi,ieegMeanSyllableShuffleRoi,1000,1,1.645);
 [zValsRawAct, pValsRaw, actClust] = timePermClusterAfterPerm(ieegMeanRoi,ieegShuffleRoi,1,1.645);

pmask = zeros(1,size(ieegMeanRoi,2));
% [~, ~, actClust]=timePermClusterZOptim(acctimeall,acctimeshuffle,10000,1,norminv(1-pVal2Cutoff));

for iClust=1:length(actClust.Size)
    if actClust.Size{iClust}>actClust.perm95
        pmask(actClust.Start{iClust}: ...
            actClust.Start{iClust}+(actClust.Size{iClust}-1)) ...
            = 1;
    end
end
visTimePlot(timeRange,ieegMeanRoiSubj,colval = colors(3,:), labels = xlab)
hold on;
visTimePlot(timeRange,squeeze((mean(ieegShuffMod(:,roiIds,:),1))),colval = [0.5 0.5 0.5], labels = xlab)
scatter(timeGammaMod(logical(pmask)),ylineval.*ones(1,sum(pmask)),'filled',MarkerEdgeColor=colors(3,:),MarkerFaceColor=colors(3,:));

groupingIds(roiIds) = 3;
title(['rMFG ' num2str(sum(roiIds))])
ieegRoiEpochScaled = ieegRoi.*pmask; [maxAcc,maxId] = max(ieegRoiEpochScaled');
timeSyllableMax{2} = timeGammaMod(maxId);
time_centroid = [];
for iElec = 1:size(ieegRoiEpochScaled,1)
    time_centroid = [time_centroid calculateTimeCentroid(timeGammaMod,ieegRoiEpochScaled(iElec,:))];
   
end
time_centroid = time_centroid(time_centroid>-0.9);
timeSyllableCentroid{2} = time_centroid;
if(showBoxPlot&&sum(pmask)>0)
    boxchart(boxplace*ones(size(time_centroid)),...
        double(time_centroid),'orientation','horizontal','BoxFaceColor',...
        colors(3,:)*0.75,'BoxWidth',boxStd, 'BoxLineColor',...
        colors(3,:),'MarkerStyle','none')
end
ylim(yrange)
pmaskroi(2,:) = pmask;

subplot(p(1),p(2),3)
bnLabels = extractBNlabels('cmfg',hemisphere)
roiIds = extractRoiIds(Subject, channelNamePooledSelect,bnLabels);
channelInfoAll = extractChannelLocation(Subject,channelNamePooledSelect);
% xyzChannel = reshape(extractfield(channelInfoAll,'xyz'),3,length(channelNamePooled));
% postAcIds = xyzChannel(2,:)<0;
% roiIds = roiIds & postAcIds;
ieegRoi = squeeze(mean(ieegMeanMod(:,roiIds,:),1));

trials2select = randperm(size(ieegMeanMod,1),groupNum);
ieegMeanRoi = squeeze(mean(mean(ieegMeanMod(:,roiIds,:),2),1))';
ieegMeanRoiTrial = squeeze(mean(ieegMeanMod(trials2select,roiIds,:),1));
[ieegMeanRoi,ieegMeanRoiSubj] = extractMeanControlSubject(ieegMeanRoiTrial,channelNamePooledSelect(roiIds));
ieegShuffleRoi = [];
for iTer = 1:1000
    trials2select = randperm(size(ieegMeanMod,1),groupNum);
    % ieegShuffleRoi(iTer, :) = squeeze(mean(mean(ieegShuffMod(trials2select,roiIds,:),2),1))';
   ieegShuffRoiTrial = squeeze(mean(ieegShuffMod(trials2select,roiIds,:),1));
    ieegShuffleRoi(iTer, :) = extractMeanControlSubject(ieegShuffRoiTrial,channelNamePooledSelect(roiIds));
end
ieegMeanSyllableRoi = squeeze(mean(ieegMeanSyllable(:,roiIds,timeSelect),2));
ieegMeanSyllableShuffleRoi = squeeze(mean(ieegMeanSyllableShuffle(:,roiIds,timeSelect),2));
% [zValsRawAct, pValsRaw, actClust] = timePermCluster(ieegMeanSyllableRoi,ieegMeanSyllableShuffleRoi,1000,1,1.645);
 [zValsRawAct, pValsRaw, actClust] = timePermClusterAfterPerm(ieegMeanRoi,ieegShuffleRoi,1,1.645);
pmask = zeros(1,size(ieegMeanRoi,2));
% [~, ~, actClust]=timePermClusterZOptim(acctimeall,acctimeshuffle,10000,1,norminv(1-pVal2Cutoff));

for iClust=1:length(actClust.Size)
    if actClust.Size{iClust}>actClust.perm95
        pmask(actClust.Start{iClust}: ...
            actClust.Start{iClust}+(actClust.Size{iClust}-1)) ...
            = 1;
    end
end
visTimePlot(timeRange,ieegMeanRoiSubj,colval = colors(10,:), labels = xlab)
hold on;
visTimePlot(timeRange,squeeze((mean(ieegShuffMod(:,roiIds,:),1))),colval = [0.5 0.5 0.5], labels = xlab)
scatter(timeGammaMod(logical(pmask)),ylineval.*ones(1,sum(pmask)),'filled',MarkerEdgeColor=colors(10,:),MarkerFaceColor=colors(10,:));

groupingIds(roiIds) = 10;
title(['cMFG ' num2str(sum(roiIds))])
ieegRoiEpochScaled = ieegRoi.*pmask; [maxAcc,maxId] = max(ieegRoiEpochScaled');
timeSyllableMax{3} = timeGammaMod(maxId);
time_centroid = [];
for iElec = 1:size(ieegRoiEpochScaled,1)
    time_centroid = [time_centroid calculateTimeCentroid(timeGammaMod,ieegRoiEpochScaled(iElec,:))];
   
end
time_centroid = time_centroid(time_centroid>-0.9);
timeSyllableCentroid{3} = time_centroid;
if(showBoxPlot&&sum(pmask)>0)
    boxchart(boxplace*ones(size(time_centroid)),...
        double(time_centroid),'orientation','horizontal','BoxFaceColor',...
        colors(10,:)*0.75,'BoxWidth',boxStd, 'BoxLineColor',...
        colors(10,:),'MarkerStyle','none')
end
ylim(yrange)
pmaskroi(3,:) = pmask;

subplot(p(1),p(2),4)
bnLabels = extractBNlabels('smc',hemisphere)
roiIds = extractRoiIds(Subject, channelNamePooledSelect,bnLabels);
ieegRoi = squeeze(mean(ieegMeanMod(:,roiIds,:),1));

trials2select = randperm(size(ieegMeanMod,1),groupNum);
ieegMeanRoi = squeeze(mean(mean(ieegMeanMod(:,roiIds,:),2),1))';
ieegMeanRoiTrial = squeeze(mean(ieegMeanMod(trials2select,roiIds,:),1));
[ieegMeanRoi,ieegMeanRoiSubj] = extractMeanControlSubject(ieegMeanRoiTrial,channelNamePooledSelect(roiIds));
ieegShuffleRoi = [];
for iTer = 1:1000
    trials2select = randperm(size(ieegMeanMod,1),groupNum);
   % ieegShuffleRoi(iTer, :) = squeeze(mean(mean(ieegShuffMod(trials2select,roiIds,:),2),1))';
   ieegShuffRoiTrial = squeeze(mean(ieegShuffMod(trials2select,roiIds,:),1));
    ieegShuffleRoi(iTer, :) = extractMeanControlSubject(ieegShuffRoiTrial,channelNamePooledSelect(roiIds));
end
ieegMeanSyllableRoi = squeeze(mean(ieegMeanSyllable(:,roiIds,timeSelect),2));
ieegMeanSyllableShuffleRoi = squeeze(mean(ieegMeanSyllableShuffle(:,roiIds,timeSelect),2));
% [zValsRawAct, pValsRaw, actClust] = timePermCluster(ieegMeanSyllableRoi,ieegMeanSyllableShuffleRoi,1000,1,1.645);
 [zValsRawAct, pValsRaw, actClust] = timePermClusterAfterPerm(ieegMeanRoi,ieegShuffleRoi,1,1.645);
pmask = zeros(1,size(ieegMeanRoi,2));

for iClust=1:length(actClust.Size)
    if actClust.Size{iClust}>actClust.perm95
        pmask(actClust.Start{iClust}: ...
            actClust.Start{iClust}+(actClust.Size{iClust}-1)) ...
            = 1;
    end
end
visTimePlot(timeRange,ieegMeanRoiSubj,colval = colors(5,:), labels = xlab)
hold on;
visTimePlot(timeRange,squeeze((mean(ieegShuffMod(:,roiIds,:),1))),colval = [0.5 0.5 0.5], labels = xlab)
scatter(timeGammaMod(logical(pmask)),ylineval.*ones(1,sum(pmask)),'filled',MarkerEdgeColor=colors(5,:),MarkerFaceColor=colors(5,:));

groupingIds(roiIds) = 5;
title(['SMC ' num2str(sum(roiIds))])
ieegRoiEpochScaled = ieegRoi.*pmask; [maxAcc,maxId] = max(ieegRoiEpochScaled');
timeSyllableMax{4} = timeGammaMod(maxId);
time_centroid = [];
for iElec = 1:size(ieegRoiEpochScaled,1)
    time_centroid = [time_centroid calculateTimeCentroid(timeGammaMod,ieegRoiEpochScaled(iElec,:))];
   
end
time_centroid = time_centroid(time_centroid>-0.9);
timeSyllableCentroid{4} = time_centroid;
if(showBoxPlot&&sum(pmask)>0)
    boxchart(boxplace*ones(size(time_centroid)),...
        double(time_centroid),'orientation','horizontal','BoxFaceColor',...
        colors(5,:)*0.75,'BoxWidth',boxStd, 'BoxLineColor',...
        colors(5,:),'MarkerStyle','none')
end
ylim(yrange)
pmaskroi(4,:) = pmask;

subplot(p(1),p(2),5)
bnLabels = extractBNlabels('ipc',hemisphere)
roiIds = extractRoiIds(Subject, channelNamePooledSelect,bnLabels);
ieegRoi = squeeze(mean(ieegMeanMod(:,roiIds,:),1));

trials2select = randperm(size(ieegMeanMod,1),groupNum);
ieegMeanRoi = squeeze(mean(mean(ieegMeanMod(:,roiIds,:),2),1))';
ieegMeanRoiTrial = squeeze(mean(ieegMeanMod(trials2select,roiIds,:),1));
[ieegMeanRoi,ieegMeanRoiSubj] = extractMeanControlSubject(ieegMeanRoiTrial,channelNamePooledSelect(roiIds));
ieegShuffleRoi = [];
for iTer = 1:1000
    trials2select = randperm(size(ieegMeanMod,1),groupNum);
    % ieegShuffleRoi(iTer, :) = squeeze(mean(mean(ieegShuffMod(trials2select,roiIds,:),2),1))';
   ieegShuffRoiTrial = squeeze(mean(ieegShuffMod(trials2select,roiIds,:),1));
    ieegShuffleRoi(iTer, :) = extractMeanControlSubject(ieegShuffRoiTrial,channelNamePooledSelect(roiIds));
end
ieegMeanSyllableRoi = squeeze(mean(ieegMeanSyllable(:,roiIds,:),2));
ieegMeanSyllableShuffleRoi = squeeze(mean(ieegMeanSyllableShuffle(:,roiIds,:),2));
% [zValsRawAct, pValsRaw, actClust] = timePermCluster(ieegMeanSyllableRoi,ieegMeanSyllableShuffleRoi,1000,1,1.645);
 [zValsRawAct, pValsRaw, actClust] = timePermClusterAfterPerm(ieegMeanRoi,ieegShuffleRoi,1,1.645);
pmask = zeros(1,size(ieegMeanRoi,2));
% [~, ~, actClust]=timePermClusterZOptim(acctimeall,acctimeshuffle,10000,1,norminv(1-pVal2Cutoff));

for iClust=1:length(actClust.Size)
    if actClust.Size{iClust}>actClust.perm95
        pmask(actClust.Start{iClust}: ...
            actClust.Start{iClust}+(actClust.Size{iClust}-1)) ...
            = 1;
    end
end
visTimePlot(timeRange,ieegMeanRoiSubj,colval = colors(2,:), labels = xlab)
hold on;
visTimePlot(timeRange,squeeze((mean(ieegShuffMod(:,roiIds,:),1))),colval = [0.5 0.5 0.5], labels = xlab)
scatter(timeGammaMod(logical(pmask)),ylineval.*ones(1,sum(pmask)),'filled',MarkerEdgeColor=colors(2,:),MarkerFaceColor=colors(2,:));

groupingIds(roiIds) = 2;
title(['IPC ' num2str(sum(roiIds))])
ieegRoiEpochScaled = ieegRoi.*pmask; [maxAcc,maxId] = max(ieegRoiEpochScaled');
timeSyllableMax{5} = timeGammaMod(maxId);
time_centroid = [];
for iElec = 1:size(ieegRoiEpochScaled,1)
    time_centroid = [time_centroid calculateTimeCentroid(timeGammaMod,ieegRoiEpochScaled(iElec,:))];
   
end
time_centroid = time_centroid(time_centroid>-0.9);
timeSyllableCentroid{5} = time_centroid;
if(showBoxPlot&&sum(pmask)>0)
    boxchart(boxplace*ones(size(time_centroid)),...
        double(time_centroid),'orientation','horizontal','BoxFaceColor',...
        colors(2,:)*0.75,'BoxWidth',boxStd, 'BoxLineColor',...
        colors(2,:),'MarkerStyle','none')
end
ylim(yrange)
pmaskroi(5,:) = pmask;

subplot(p(1),p(2),6)
bnLabels = extractBNlabels('insula',hemisphere)
roiIds = extractRoiIds(Subject, channelNamePooledSelect,bnLabels);
ieegRoi = squeeze(mean(ieegMeanMod(:,roiIds,:),1));

trials2select = randperm(size(ieegMeanMod,1),groupNum);
ieegMeanRoi = squeeze(mean(mean(ieegMeanMod(:,roiIds,:),2),1))';
ieegMeanRoiTrial = squeeze(mean(ieegMeanMod(trials2select,roiIds,:),1));
[ieegMeanRoi,ieegMeanRoiSubj] = extractMeanControlSubject(ieegMeanRoiTrial,channelNamePooledSelect(roiIds));
ieegShuffleRoi = [];
for iTer = 1:1000
    trials2select = randperm(size(ieegMeanMod,1),groupNum);
    % ieegShuffleRoi(iTer, :) = squeeze(mean(mean(ieegShuffMod(trials2select,roiIds,:),2),1))';
   ieegShuffRoiTrial = squeeze(mean(ieegShuffMod(trials2select,roiIds,:),1));
    ieegShuffleRoi(iTer, :) = extractMeanControlSubject(ieegShuffRoiTrial,channelNamePooledSelect(roiIds));
end
ieegMeanSyllableRoi = squeeze(mean(ieegMeanSyllable(:,roiIds,:),2));
ieegMeanSyllableShuffleRoi = squeeze(mean(ieegMeanSyllableShuffle(:,roiIds,:),2));
%[zValsRawAct, pValsRaw, actClust] = timePermCluster(ieegMeanSyllableRoi,ieegMeanSyllableShuffleRoi,1000,1,1.645);
[zValsRawAct, pValsRaw, actClust] = timePermClusterAfterPerm(ieegMeanRoi,ieegShuffleRoi,1,1.645);
pmask = zeros(1,size(ieegMeanRoi,2));
% [~, ~, actClust]=timePermClusterZOptim(acctimeall,acctimeshuffle,10000,1,norminv(1-pVal2Cutoff));

for iClust=1:length(actClust.Size)
    if actClust.Size{iClust}>actClust.perm95
        pmask(actClust.Start{iClust}: ...
            actClust.Start{iClust}+(actClust.Size{iClust}-1)) ...
            = 1;
    end
end
visTimePlot(timeRange,ieegMeanRoiSubj,colval = colors(6,:), labels = xlab)
hold on;
visTimePlot(timeRange,squeeze((mean(ieegShuffMod(:,roiIds,:),1))),colval = [0.5 0.5 0.5], labels = xlab)
scatter(timeGammaMod(logical(pmask)),ylineval.*ones(1,sum(pmask)),'filled',MarkerEdgeColor=colors(6,:),MarkerFaceColor=colors(6,:));

groupingIds(roiIds) = 1;
title(['Insula ' num2str(sum(roiIds))])
ieegRoiEpochScaled = ieegRoi.*pmask; [maxAcc,maxId] = max(ieegRoiEpochScaled');
timeSyllableMax{6} = timeGammaMod(maxId);
time_centroid = [];
for iElec = 1:size(ieegRoiEpochScaled,1)
    time_centroid = [time_centroid calculateTimeCentroid(timeGammaMod,ieegRoiEpochScaled(iElec,:))];
   
end
time_centroid = time_centroid(time_centroid>-0.9);
timeSyllableCentroid{6} = time_centroid;
if(showBoxPlot&&sum(pmask)>0)
    boxchart(boxplace*ones(size(time_centroid)),...
        double(time_centroid),'orientation','horizontal','BoxFaceColor',...
        colors(6,:)*0.75,'BoxWidth',boxStd, 'BoxLineColor',...
        colors(6,:),'MarkerStyle','none')
end
ylim(yrange)
pmaskroi(6,:) = pmask;


subplot(p(1),p(2),7)
bnLabels = extractBNlabels('astg',hemisphere)
roiIds = extractRoiIds(Subject, channelNamePooledSelect,bnLabels);
ieegRoi = squeeze(mean(ieegMeanMod(:,roiIds,:),1));

trials2select = randperm(size(ieegMeanMod,1),groupNum);
ieegMeanRoi = squeeze(mean(mean(ieegMeanMod(:,roiIds,:),2),1))';
ieegMeanRoiTrial = squeeze(mean(ieegMeanMod(trials2select,roiIds,:),1));
[ieegMeanRoi,ieegMeanRoiSubj] = extractMeanControlSubject(ieegMeanRoiTrial,channelNamePooledSelect(roiIds));
ieegShuffleRoi = [];
for iTer = 1:1000
    trials2select = randperm(size(ieegMeanMod,1),groupNum);
   % ieegShuffleRoi(iTer, :) = squeeze(mean(mean(ieegShuffMod(trials2select,roiIds,:),2),1))';
   ieegShuffRoiTrial = squeeze(mean(ieegShuffMod(trials2select,roiIds,:),1));
    ieegShuffleRoi(iTer, :) = extractMeanControlSubject(ieegShuffRoiTrial,channelNamePooledSelect(roiIds));
end
ieegMeanSyllableRoi = squeeze(mean(ieegMeanSyllable(:,roiIds,:),2));
ieegMeanSyllableShuffleRoi = squeeze(mean(ieegMeanSyllableShuffle(:,roiIds,:),2));
% [zValsRawAct, pValsRaw, actClust] = timePermCluster(ieegMeanSyllableRoi,ieegMeanSyllableShuffleRoi,1000,1,1.645);
 [zValsRawAct, pValsRaw, actClust] = timePermClusterAfterPerm(ieegMeanRoi,ieegShuffleRoi,1,1.645);
pmask = zeros(1,size(ieegMeanRoi,2));
% [~, ~, actClust]=timePermClusterZOptim(acctimeall,acctimeshuffle,10000,1,norminv(1-pVal2Cutoff));

for iClust=1:length(actClust.Size)
    if actClust.Size{iClust}>actClust.perm95
        pmask(actClust.Start{iClust}: ...
            actClust.Start{iClust}+(actClust.Size{iClust}-1)) ...
            = 1;
    end
end
visTimePlot(timeRange,ieegMeanRoiSubj,colval = colors(9,:), labels = xlab)
hold on;
visTimePlot(timeRange,squeeze((mean(ieegShuffMod(:,roiIds,:),1))),colval = [0.5 0.5 0.5], labels = xlab)
scatter(timeGammaMod(logical(pmask)),ylineval.*ones(1,sum(pmask)),'filled',MarkerEdgeColor=colors(9,:),MarkerFaceColor=colors(9,:));

groupingIds(roiIds) = 9;
title(['aSTG ' num2str(sum(roiIds))])
ieegRoiEpochScaled = ieegRoi.*pmask; [maxAcc,maxId] = max(ieegRoiEpochScaled');
timeSyllableMax{7} = timeGammaMod(maxId);
time_centroid = [];
for iElec = 1:size(ieegRoiEpochScaled,1)
    time_centroid = [time_centroid calculateTimeCentroid(timeGammaMod,ieegRoiEpochScaled(iElec,:))];
   
end
time_centroid = time_centroid(time_centroid>-0.9);
timeSyllableCentroid{7} = time_centroid;
if(showBoxPlot&&sum(pmask)>0)
    boxchart(boxplace*ones(size(time_centroid)),...
        double(time_centroid),'orientation','horizontal','BoxFaceColor',...
        colors(9,:)*0.75,'BoxWidth',boxStd, 'BoxLineColor',...
        colors(9,:),'MarkerStyle','none')
end
ylim(yrange)
pmaskroi(7,:) = pmask;
subplot(p(1),p(2),8)
bnLabels = extractBNlabels('pstg',hemisphere)
roiIds = extractRoiIds(Subject, channelNamePooledSelect,bnLabels);
ieegRoi = squeeze(mean(ieegMeanMod(:,roiIds,:),1));

trials2select = randperm(size(ieegMeanMod,1),groupNum);
ieegMeanRoi = squeeze(mean(mean(ieegMeanMod(:,roiIds,:),2),1))';
ieegMeanRoiTrial = squeeze(mean(ieegMeanMod(trials2select,roiIds,:),1));
[ieegMeanRoi,ieegMeanRoiSubj] = extractMeanControlSubject(ieegMeanRoiTrial,channelNamePooledSelect(roiIds));
ieegShuffleRoi = [];
for iTer = 1:1000
    trials2select = randperm(size(ieegMeanMod,1),groupNum);
    % ieegShuffleRoi(iTer, :) = squeeze(mean(mean(ieegShuffMod(trials2select,roiIds,:),2),1))';
   ieegShuffRoiTrial = squeeze(mean(ieegShuffMod(trials2select,roiIds,:),1));
    ieegShuffleRoi(iTer, :) = extractMeanControlSubject(ieegShuffRoiTrial,channelNamePooledSelect(roiIds));
end
ieegMeanSyllableRoi = squeeze(mean(ieegMeanSyllable(:,roiIds,:),2));
ieegMeanSyllableShuffleRoi = squeeze(mean(ieegMeanSyllableShuffle(:,roiIds,:),2));
% [zValsRawAct, pValsRaw, actClust] = timePermCluster(ieegMeanSyllableRoi,ieegMeanSyllableShuffleRoi,1000,1,1.645);
 [zValsRawAct, pValsRaw, actClust] = timePermClusterAfterPerm(ieegMeanRoi,ieegShuffleRoi,1,1.645);
pmask = zeros(1,size(ieegMeanRoi,2));
% [~, ~, actClust]=timePermClusterZOptim(acctimeall,acctimeshuffle,10000,1,norminv(1-pVal2Cutoff));

for iClust=1:length(actClust.Size)
    if actClust.Size{iClust}>actClust.perm95
        pmask(actClust.Start{iClust}: ...
            actClust.Start{iClust}+(actClust.Size{iClust}-1)) ...
            = 1;
    end
end
visTimePlot(timeRange,ieegMeanRoiSubj,colval = colors(1,:), labels = xlab)
hold on;
visTimePlot(timeRange,squeeze((mean(ieegShuffMod(:,roiIds,:),1))),colval = [0.5 0.5 0.5], labels = xlab)
scatter(timeGammaMod(logical(pmask)),ylineval.*ones(1,sum(pmask)),'filled',MarkerEdgeColor=colors(1,:),MarkerFaceColor=colors(1,:));

groupingIds(roiIds) = 1;
ieegRoiEpochScaled = ieegRoi.*pmask; [maxAcc,maxId] = max(ieegRoiEpochScaled');
timeSyllableMax{8} = timeGammaMod(maxId);
time_centroid = [];
for iElec = 1:size(ieegRoiEpochScaled,1)
    time_centroid = [time_centroid calculateTimeCentroid(timeGammaMod,ieegRoiEpochScaled(iElec,:))];
   
end
time_centroid = time_centroid(time_centroid>-0.9);
timeSyllableCentroid{8} = time_centroid;
if(showBoxPlot&&sum(pmask)>0)
    boxchart(boxplace*ones(size(time_centroid)),...
        double(time_centroid),'orientation','horizontal','BoxFaceColor',...
        colors(1,:)*0.75,'BoxWidth',boxStd, 'BoxLineColor',...
        colors(1,:),'MarkerStyle','none')
end
title(['pSTG ' num2str(sum(roiIds))])
ylim(yrange)
pmaskroi(8,:) = pmask;
subplot(p(1),p(2),9)
bnLabels = extractBNlabels('sts',hemisphere)
roiIds = extractRoiIds(Subject, channelNamePooledSelect,bnLabels);
ieegRoi = squeeze(mean(ieegMeanMod(:,roiIds,:),1));

trials2select = randperm(size(ieegMeanMod,1),groupNum);
ieegMeanRoi = squeeze(mean(mean(ieegMeanMod(:,roiIds,:),2),1))';
ieegMeanRoiTrial = squeeze(mean(ieegMeanMod(trials2select,roiIds,:),1));
[ieegMeanRoi,ieegMeanRoiSubj] = extractMeanControlSubject(ieegMeanRoiTrial,channelNamePooledSelect(roiIds));
ieegShuffleRoi = [];
for iTer = 1:1000
    trials2select = randperm(size(ieegMeanMod,1),groupNum);
    % ieegShuffleRoi(iTer, :) = squeeze(mean(mean(ieegShuffMod(trials2select,roiIds,:),2),1))';
   ieegShuffRoiTrial = squeeze(mean(ieegShuffMod(trials2select,roiIds,:),1));
    ieegShuffleRoi(iTer, :) = extractMeanControlSubject(ieegShuffRoiTrial,channelNamePooledSelect(roiIds));
end
ieegMeanSyllableRoi = squeeze(mean(ieegMeanSyllable(:,roiIds,:),2));
ieegMeanSyllableShuffleRoi = squeeze(mean(ieegMeanSyllableShuffle(:,roiIds,:),2));
% [zValsRawAct, pValsRaw, actClust] = timePermCluster(ieegMeanSyllableRoi,ieegMeanSyllableShuffleRoi,1000,1,1.645);
 [zValsRawAct, pValsRaw, actClust] = timePermClusterAfterPerm(ieegMeanRoi,ieegShuffleRoi,1,1.645);
pmask = zeros(1,size(ieegMeanRoi,2));
% [~, ~, actClust]=timePermClusterZOptim(acctimeall,acctimeshuffle,10000,1,norminv(1-pVal2Cutoff));

for iClust=1:length(actClust.Size)
    if actClust.Size{iClust}>actClust.perm95
        pmask(actClust.Start{iClust}: ...
            actClust.Start{iClust}+(actClust.Size{iClust}-1)) ...
            = 1;
    end
end
visTimePlot(timeRange,ieegMeanRoiSubj,colval = colors(8,:), labels = xlab)
hold on;
visTimePlot(timeRange,squeeze((mean(ieegShuffMod(:,roiIds,:),1))),colval = [0.5 0.5 0.5], labels = xlab)
scatter(timeGammaMod(logical(pmask)),ylineval.*ones(1,sum(pmask)),'filled',MarkerEdgeColor=colors(8,:),MarkerFaceColor=colors(8,:));

groupingIds(roiIds) = 1;
ieegRoiEpochScaled = ieegRoi.*pmask; [maxAcc,maxId] = max(ieegRoiEpochScaled');
timeSyllableMax{9} = timeGammaMod(maxId);
time_centroid = [];
for iElec = 1:size(ieegRoiEpochScaled,1)
    time_centroid = [time_centroid calculateTimeCentroid(timeGammaMod,ieegRoiEpochScaled(iElec,:))];
   
end
time_centroid = time_centroid(time_centroid>-0.9);
timeSyllableCentroid{9} = time_centroid;
if(showBoxPlot&&sum(pmask)>0)
    boxchart(boxplace*ones(size(time_centroid)),...
        double(time_centroid),'orientation','horizontal','BoxFaceColor',...
        colors(8,:)*0.75,'BoxWidth',boxStd, 'BoxLineColor',...
        colors(8,:),'MarkerStyle','none')
end
title(['STS ' num2str(sum(roiIds))])
ylim(yrange)
pmaskroi(9,:) = pmask;

subplot(p(1),p(2),10)
bnLabels = extractBNlabels('heschl',hemisphere)

roiIds = extractRoiIds(Subject, channelNamePooledSelect,bnLabels);
% xyzChannel = reshape(extractfield(channelInfoAll,'xyz'),3,length(channelNamePooled));
% postAcIds = xyzChannel(2,:)>0;
%roiIds = roiIds & preAcIds;

ieegRoi = squeeze(mean(ieegMeanMod(:,roiIds,:),1));

trials2select = randperm(size(ieegMeanMod,1),groupNum);
ieegMeanRoi = squeeze(mean(mean(ieegMeanMod(:,roiIds,:),2),1))';
ieegMeanRoiTrial = squeeze(mean(ieegMeanMod(trials2select,roiIds,:),1));
[ieegMeanRoi,ieegMeanRoiSubj] = extractMeanControlSubject(ieegMeanRoiTrial,channelNamePooledSelect(roiIds));
ieegShuffleRoi = [];
for iTer = 1:1000
    trials2select = randperm(size(ieegMeanMod,1),groupNum);
    % ieegShuffleRoi(iTer, :) = squeeze(mean(mean(ieegShuffMod(trials2select,roiIds,:),2),1))';
   ieegShuffRoiTrial = squeeze(mean(ieegShuffMod(trials2select,roiIds,:),1));
    ieegShuffleRoi(iTer, :) = extractMeanControlSubject(ieegShuffRoiTrial,channelNamePooledSelect(roiIds));
end
ieegMeanSyllableRoi = squeeze(mean(ieegMeanSyllable(:,roiIds,timeSelect),2));
ieegMeanSyllableShuffleRoi = squeeze(mean(ieegMeanSyllableShuffle(:,roiIds,timeSelect),2));
% [zValsRawAct, pValsRaw, actClust] = timePermCluster(ieegMeanSyllableRoi,ieegMeanSyllableShuffleRoi,1000,1,1.645);
 [zValsRawAct, pValsRaw, actClust] = timePermClusterAfterPerm(ieegMeanRoi,ieegShuffleRoi,1,1.645);
pmask = zeros(1,size(ieegMeanRoi,2));
% [~, ~, actClust]=timePermClusterZOptim(acctimeall,acctimeshuffle,10000,1,norminv(1-pVal2Cutoff));

for iClust=1:length(actClust.Size)
    if actClust.Size{iClust}>actClust.perm95
        pmask(actClust.Start{iClust}: ...
            actClust.Start{iClust}+(actClust.Size{iClust}-1)) ...
            = 1;
    end
end
visTimePlot(timeRange,ieegMeanRoiSubj,colval = colors(7,:), labels = xlab )
hold on;
visTimePlot(timeRange,squeeze((mean(ieegShuffMod(:,roiIds,:),1))),colval = [0.5 0.5 0.5], labels = xlab)
scatter(timeGammaMod(logical(pmask)),ylineval.*ones(1,sum(pmask)),'filled',MarkerEdgeColor=colors(7,:),MarkerFaceColor=colors(7,:));
groupingIds(roiIds) = 7;
title(['Heschl ' num2str(sum(roiIds))])
ieegRoiEpochScaled = ieegRoi.*pmask; [maxAcc,maxId] = max(ieegRoiEpochScaled');

timeSyllableMax{10} = timeGammaMod(maxId);
time_centroid = [];
for iElec = 1:size(ieegRoiEpochScaled,1)
    time_centroid = [time_centroid calculateTimeCentroid(timeGammaMod,ieegRoiEpochScaled(iElec,:))];
   
end
time_centroid = time_centroid(time_centroid>-0.9);
timeSyllableCentroid{10} = time_centroid;
if(showBoxPlot&&sum(pmask)>0)
    boxchart(boxplace*ones(size(time_centroid)),...
        double(time_centroid),'orientation','horizontal','BoxFaceColor',...
        colors(7,:)*0.75,'BoxWidth',boxStd, 'BoxLineColor',...
        colors(7,:),'MarkerStyle','none')
end
ylim(yrange)
pmaskroi(10,:) = pmask;


% visTimePlot(timeRange,modelweightnorm(groupingIds==0,:),colval = colors(9,:))
%  groupingIds(groupingIds==0) = 9;

% xlabel('Time from response onset')
% ylabel('Beta')
sgtitle('')
%%
load colors.mat



hemisphere = 'lh';

figure;
bnLabels = extractBNlabels('ifg',hemisphere)
roiIds = extractRoiIds(Subject, channelNamePooledSelect,bnLabels);
ieegRoi = squeeze(mean(ieegMeanMod(:,roiIds,:),1));

trials2select = randperm(size(ieegMeanMod,1),groupNum);

ieegMeanRoiTrial = squeeze(mean(ieegMeanMod(trials2select,roiIds,:),1));
[ieegMeanRoi,ieegMeanRoiSubj] = extractMeanControlSubject(ieegMeanRoiTrial,channelNamePooledSelect(roiIds));

visTimePlot(timeRange,ieegMeanRoiSubj,colval = colors(4,:), labels = xlab)
scatter(timeGammaMod(logical(pmaskroi(1,:))),ylineval.*ones(1,sum(pmaskroi(1,:))),'filled',MarkerEdgeColor=colors(4,:),MarkerFaceColor=colors(4,:));
ylim(yrange)


hold on;
bnLabels = extractBNlabels('rmfg',hemisphere)
roiIds = extractRoiIds(Subject, channelNamePooledSelect,bnLabels);


trials2select = randperm(size(ieegMeanMod,1),groupNum);

ieegMeanRoiTrial = squeeze(mean(ieegMeanMod(trials2select,roiIds,:),1));
[ieegMeanRoi,ieegMeanRoiSubj] = extractMeanControlSubject(ieegMeanRoiTrial,channelNamePooledSelect(roiIds));

visTimePlot(timeRange,ieegMeanRoiSubj,colval = colors(3,:), labels = xlab)
scatter(timeGammaMod(logical(pmaskroi(2,:))),ylineval.*ones(1,sum(pmaskroi(2,:))),'filled',MarkerEdgeColor=colors(3,:),MarkerFaceColor=colors(3,:));


hold on;
bnLabels = extractBNlabels('cmfg',hemisphere)
roiIds = extractRoiIds(Subject, channelNamePooledSelect,bnLabels);


trials2select = randperm(size(ieegMeanMod,1),groupNum);
ieegMeanRoi = squeeze(mean(mean(ieegMeanMod(:,roiIds,:),2),1))';
ieegMeanRoiTrial = squeeze(mean(ieegMeanMod(trials2select,roiIds,:),1));
[ieegMeanRoi,ieegMeanRoiSubj] = extractMeanControlSubject(ieegMeanRoiTrial,channelNamePooledSelect(roiIds));

visTimePlot(timeRange,ieegMeanRoiSubj,colval = colors(10,:), labels = xlab)
hold on;
scatter(timeGammaMod(logical(pmaskroi(3,:))),ylineval.*ones(1,sum(pmaskroi(3,:))),'filled',MarkerEdgeColor=colors(10,:),MarkerFaceColor=colors(10,:));
ylim(yrange)


figure;

bnLabels = extractBNlabels('smc',hemisphere)
roiIds = extractRoiIds(Subject, channelNamePooledSelect,bnLabels);
ieegRoi = squeeze(mean(ieegMeanMod(:,roiIds,:),1));

trials2select = randperm(size(ieegMeanMod,1),groupNum);
ieegMeanRoi = squeeze(mean(mean(ieegMeanMod(:,roiIds,:),2),1))';
ieegMeanRoiTrial = squeeze(mean(ieegMeanMod(trials2select,roiIds,:),1));
[ieegMeanRoi,ieegMeanRoiSubj] = extractMeanControlSubject(ieegMeanRoiTrial,channelNamePooledSelect(roiIds));

visTimePlot(timeRange,ieegMeanRoiSubj,colval = colors(5,:), labels = xlab)
hold on;
scatter(timeGammaMod(logical(pmaskroi(4,:))),ylineval.*ones(1,sum(pmaskroi(4,:))),'filled',MarkerEdgeColor=colors(5,:),MarkerFaceColor=colors(5,:));
ylim(yrange)


hold on;
bnLabels = extractBNlabels('ipc',hemisphere)
roiIds = extractRoiIds(Subject, channelNamePooledSelect,bnLabels);
ieegRoi = squeeze(mean(ieegMeanMod(:,roiIds,:),1));

trials2select = randperm(size(ieegMeanMod,1),groupNum);
ieegMeanRoi = squeeze(mean(mean(ieegMeanMod(:,roiIds,:),2),1))';
ieegMeanRoiTrial = squeeze(mean(ieegMeanMod(trials2select,roiIds,:),1));
[ieegMeanRoi,ieegMeanRoiSubj] = extractMeanControlSubject(ieegMeanRoiTrial,channelNamePooledSelect(roiIds));

visTimePlot(timeRange,ieegMeanRoiSubj,colval = colors(2,:), labels = xlab)
scatter(timeGammaMod(logical(pmaskroi(5,:))),ylineval.*ones(1,sum(pmaskroi(5,:))),'filled',MarkerEdgeColor=colors(2,:),MarkerFaceColor=colors(2,:));
ylim(yrange)


figure;
bnLabels = extractBNlabels('insula',hemisphere)
roiIds = extractRoiIds(Subject, channelNamePooledSelect,bnLabels);
ieegRoi = squeeze(mean(ieegMeanMod(:,roiIds,:),1));

trials2select = randperm(size(ieegMeanMod,1),groupNum);
ieegMeanRoi = squeeze(mean(mean(ieegMeanMod(:,roiIds,:),2),1))';
ieegMeanRoiTrial = squeeze(mean(ieegMeanMod(trials2select,roiIds,:),1));
[ieegMeanRoi,ieegMeanRoiSubj] = extractMeanControlSubject(ieegMeanRoiTrial,channelNamePooledSelect(roiIds));

visTimePlot(timeRange,ieegMeanRoiSubj,colval = colors(6,:), labels = xlab)
hold on;
scatter(timeGammaMod(logical(pmaskroi(6,:))),ylineval.*ones(1,sum(pmaskroi(6,:))),'filled',MarkerEdgeColor=colors(6,:),MarkerFaceColor=colors(6,:));

ylim(yrange)



hold on;
bnLabels = extractBNlabels('astg',hemisphere)
roiIds = extractRoiIds(Subject, channelNamePooledSelect,bnLabels);
ieegRoi = squeeze(mean(ieegMeanMod(:,roiIds,:),1));

trials2select = randperm(size(ieegMeanMod,1),groupNum);
ieegMeanRoi = squeeze(mean(mean(ieegMeanMod(:,roiIds,:),2),1))';
ieegMeanRoiTrial = squeeze(mean(ieegMeanMod(trials2select,roiIds,:),1));
[ieegMeanRoi,ieegMeanRoiSubj] = extractMeanControlSubject(ieegMeanRoiTrial,channelNamePooledSelect(roiIds));

visTimePlot(timeRange,ieegMeanRoiSubj,colval = colors(9,:), labels = xlab)
hold on;
scatter(timeGammaMod(logical(pmaskroi(7,:))),ylineval.*ones(1,sum(pmaskroi(7,:))),'filled',MarkerEdgeColor=colors(9,:),MarkerFaceColor=colors(9,:));


ylim(yrange)


hold on;
bnLabels = extractBNlabels('pstg',hemisphere)
roiIds = extractRoiIds(Subject, channelNamePooledSelect,bnLabels);
ieegRoi = squeeze(mean(ieegMeanMod(:,roiIds,:),1));

trials2select = randperm(size(ieegMeanMod,1),groupNum);
ieegMeanRoi = squeeze(mean(mean(ieegMeanMod(:,roiIds,:),2),1))';
ieegMeanRoiTrial = squeeze(mean(ieegMeanMod(trials2select,roiIds,:),1));
[ieegMeanRoi,ieegMeanRoiSubj] = extractMeanControlSubject(ieegMeanRoiTrial,channelNamePooledSelect(roiIds));

visTimePlot(timeRange,ieegMeanRoiSubj,colval = colors(1,:), labels = xlab)
hold on;
scatter(timeGammaMod(logical(pmaskroi(8,:))),ylineval.*ones(1,sum(pmaskroi(8,:))),'filled',MarkerEdgeColor=colors(1,:),MarkerFaceColor=colors(1,:));
ylim(yrange)

hold on;
bnLabels = extractBNlabels('sts',hemisphere)
roiIds = extractRoiIds(Subject, channelNamePooledSelect,bnLabels);
ieegRoi = squeeze(mean(ieegMeanMod(:,roiIds,:),1));

trials2select = randperm(size(ieegMeanMod,1),groupNum);
ieegMeanRoi = squeeze(mean(mean(ieegMeanMod(:,roiIds,:),2),1))';
ieegMeanRoiTrial = squeeze(mean(ieegMeanMod(trials2select,roiIds,:),1));
[ieegMeanRoi,ieegMeanRoiSubj] = extractMeanControlSubject(ieegMeanRoiTrial,channelNamePooledSelect(roiIds));

visTimePlot(timeRange,ieegMeanRoiSubj,colval = colors(8,:), labels = xlab)
hold on;
scatter(timeGammaMod(logical(pmaskroi(9,:))),ylineval.*ones(1,sum(pmaskroi(9,:))),'filled',MarkerEdgeColor=colors(8,:),MarkerFaceColor=colors(8,:));
ylim(yrange)


hold on;
bnLabels = extractBNlabels('heschl',hemisphere)

roiIds = extractRoiIds(Subject, channelNamePooledSelect,bnLabels);
ieegRoi = squeeze(mean(ieegMeanMod(:,roiIds,:),1));

trials2select = randperm(size(ieegMeanMod,1),groupNum);
ieegMeanRoi = squeeze(mean(mean(ieegMeanMod(:,roiIds,:),2),1))';
ieegMeanRoiTrial = squeeze(mean(ieegMeanMod(trials2select,roiIds,:),1));
[ieegMeanRoi,ieegMeanRoiSubj] = extractMeanControlSubject(ieegMeanRoiTrial,channelNamePooledSelect(roiIds));

visTimePlot(timeRange,ieegMeanRoiSubj,colval = colors(7,:), labels = xlab )
hold on;
scatter(timeGammaMod(logical(pmaskroi(10,:))),ylineval.*ones(1,sum(pmaskroi(10,:))),'filled',MarkerEdgeColor=colors(7,:),MarkerFaceColor=colors(7,:));

ylim(yrange)

%%
colGroup = lines(3);
figure;


hold on;
timeStart = -0.9;

hemisphere = 'rh';

bnLabels = [extractBNlabels('ifg',hemisphere) extractBNlabels('rmfg',hemisphere) extractBNlabels('cmfg',hemisphere) ]
roiIds = extractRoiIds(Subject, channelNamePooledSelect,bnLabels);
ieegRoi = squeeze(mean(ieegMeanMod(:,roiIds,:),1));

trials2select = randperm(size(ieegMeanMod,1),groupNum);
ieegMeanRoi = squeeze(mean(mean(ieegMeanMod(:,roiIds,:),2),1))';
ieegMeanRoiTrial = squeeze(mean(ieegMeanMod(trials2select,roiIds,:),1));
[ieegMeanRoi,ieegMeanRoiSubj] = extractMeanControlSubject(ieegMeanRoiTrial,channelNamePooledSelect(roiIds));
ieegShuffleRoi = [];
for iTer = 1:1000
    trials2select = randperm(size(ieegMeanMod,1),groupNum);
   % ieegShuffleRoi(iTer, :) = squeeze(mean(mean(ieegShuffMod(trials2select,roiIds,:),2),1))';
   ieegShuffRoiTrial = squeeze(mean(ieegShuffMod(trials2select,roiIds,:),1));
    ieegShuffleRoi(iTer, :) = extractMeanControlSubject(ieegShuffRoiTrial,channelNamePooledSelect(roiIds));
end
ieegMeanSyllableRoi = squeeze(mean(ieegMeanSyllable(:,roiIds,timeSelect),2));
ieegMeanSyllableShuffleRoi = squeeze(mean(ieegMeanSyllableShuffle(:,roiIds,timeSelect),2));
%[zValsRawAct, pValsRaw, actClust] = timePermCluster(ieegMeanSyllableRoi,ieegMeanSyllableShuffleRoi,1000,1,1.645);
 [zValsRawAct, pValsRaw, actClust] = timePermClusterAfterPerm(ieegMeanRoi,ieegShuffleRoi,1,1.645);
pmask = zeros(1,size(ieegMeanRoi,2));
% [~, ~, actClust]=timePermClusterZOptim(acctimeall,acctimeshuffle,10000,1,norminv(1-pVal2Cutoff));

for iClust=1:length(actClust.Size)
    if actClust.Size{iClust}>actClust.perm95
        pmask(actClust.Start{iClust}: ...
            actClust.Start{iClust}+(actClust.Size{iClust}-1)) ...
            = 1;
    end
end
visTimePlot(timeRange,ieegMeanRoiSubj,colval = colGroup(1,:), labels = xlab)
hold on;
visTimePlot(timeRange,squeeze((mean(ieegShuffMod(:,roiIds,:),1))),colval = [0.5 0.5 0.5], labels = xlab)
scatter(timeGammaMod(logical(pmask)),ylineval.*ones(1,sum(pmask)),'filled',MarkerEdgeColor=colGroup(1,:),MarkerFaceColor=colGroup(1,:));

ieegRoiEpochScaled = ieegRoi.*pmask; [maxAcc,maxId] = max(ieegRoiEpochScaled');
timeSyllableMaxGroup{1} = timeGammaMod(maxId);
time_centroid = [];
for iElec = 1:size(ieegRoiEpochScaled,1)
    time_centroid = [time_centroid calculateTimeCentroid(timeGammaMod,ieegRoiEpochScaled(iElec,:))];
   
end
time_centroid = time_centroid(time_centroid>timeStart);
timeSyllableCentroidGroup{1} = time_centroid;
ylim(yrange)


figure;
bnLabels = [extractBNlabels('smc',hemisphere) extractBNlabels('ipc',hemisphere)]
roiIds = extractRoiIds(Subject, channelNamePooledSelect,bnLabels);
ieegRoi = squeeze(mean(ieegMeanMod(:,roiIds,:),1));

trials2select = randperm(size(ieegMeanMod,1),groupNum);
ieegMeanRoi = squeeze(mean(mean(ieegMeanMod(:,roiIds,:),2),1))';
ieegMeanRoiTrial = squeeze(mean(ieegMeanMod(trials2select,roiIds,:),1));
[ieegMeanRoi,ieegMeanRoiSubj] = extractMeanControlSubject(ieegMeanRoiTrial,channelNamePooledSelect(roiIds));
ieegShuffleRoi = [];
for iTer = 1:100
    trials2select = randperm(size(ieegMeanMod,1),groupNum);
    % ieegShuffleRoi(iTer, :) = squeeze(mean(mean(ieegShuffMod(trials2select,roiIds,:),2),1))';
   ieegShuffRoiTrial = squeeze(mean(ieegShuffMod(trials2select,roiIds,:),1));
    ieegShuffleRoi(iTer, :) = extractMeanControlSubject(ieegShuffRoiTrial,channelNamePooledSelect(roiIds));
end
ieegMeanSyllableRoi = squeeze(mean(ieegMeanSyllable(:,roiIds,timeSelect),2));
ieegMeanSyllableShuffleRoi = squeeze(mean(ieegMeanSyllableShuffle(:,roiIds,timeSelect),2));
% [zValsRawAct, pValsRaw, actClust] = timePermCluster(ieegMeanSyllableRoi,ieegMeanSyllableShuffleRoi,1000,1,1.645);
 [zValsRawAct, pValsRaw, actClust] = timePermClusterAfterPerm(ieegMeanRoi,ieegShuffleRoi,1,1.645);

pmask = zeros(1,size(ieegMeanRoi,2));
% [~, ~, actClust]=timePermClusterZOptim(acctimeall,acctimeshuffle,10000,1,norminv(1-pVal2Cutoff));

for iClust=1:length(actClust.Size)
    if actClust.Size{iClust}>actClust.perm95
        pmask(actClust.Start{iClust}: ...
            actClust.Start{iClust}+(actClust.Size{iClust}-1)) ...
            = 1;
    end
end
visTimePlot(timeRange,ieegMeanRoiSubj,colval = colGroup(2,:), labels = xlab)
hold on;
visTimePlot(timeRange,squeeze((mean(ieegShuffMod(:,roiIds,:),1))),colval = [0.5 0.5 0.5], labels = xlab)
scatter(timeGammaMod(logical(pmask)),ylineval.*ones(1,sum(pmask)),'filled',MarkerEdgeColor=colGroup(2,:),MarkerFaceColor=colGroup(2,:));

ylim(yrange)

ieegRoiEpochScaled = ieegRoi.*pmask; [maxAcc,maxId] = max(ieegRoiEpochScaled');
timeSyllableMaxGroup{2} = timeGammaMod(maxId);
time_centroid = [];
for iElec = 1:size(ieegRoiEpochScaled,1)
    time_centroid = [time_centroid calculateTimeCentroid(timeGammaMod,ieegRoiEpochScaled(iElec,:))];
   
end
time_centroid = time_centroid(time_centroid>timeStart);
timeSyllableCentroidGroup{2} = time_centroid;


figure;
bnLabels = [extractBNlabels('astg',hemisphere) extractBNlabels('pstg',hemisphere) extractBNlabels('insula',hemisphere) extractBNlabels('heschl',hemisphere) extractBNlabels('sts',hemisphere)]
roiIds = extractRoiIds(Subject, channelNamePooledSelect,bnLabels);
channelInfoAll = extractChannelLocation(Subject,channelNamePooledSelect);
% xyzChannel = reshape(extractfield(channelInfoAll,'xyz'),3,length(channelNamePooled));
% postAcIds = xyzChannel(2,:)<0;
% roiIds = roiIds & postAcIds;
ieegRoi = squeeze(mean(ieegMeanMod(:,roiIds,:),1));

trials2select = randperm(size(ieegMeanMod,1),groupNum);
ieegMeanRoi = squeeze(mean(mean(ieegMeanMod(:,roiIds,:),2),1))';
ieegMeanRoiTrial = squeeze(mean(ieegMeanMod(trials2select,roiIds,:),1));
[ieegMeanRoi,ieegMeanRoiSubj] = extractMeanControlSubject(ieegMeanRoiTrial,channelNamePooledSelect(roiIds));
ieegShuffleRoi = [];
for iTer = 1:1000
    trials2select = randperm(size(ieegMeanMod,1),groupNum);
    % ieegShuffleRoi(iTer, :) = squeeze(mean(mean(ieegShuffMod(trials2select,roiIds,:),2),1))';
   ieegShuffRoiTrial = squeeze(mean(ieegShuffMod(trials2select,roiIds,:),1));
    ieegShuffleRoi(iTer, :) = extractMeanControlSubject(ieegShuffRoiTrial,channelNamePooledSelect(roiIds));
end
ieegMeanSyllableRoi = squeeze(mean(ieegMeanSyllable(:,roiIds,timeSelect),2));
ieegMeanSyllableShuffleRoi = squeeze(mean(ieegMeanSyllableShuffle(:,roiIds,timeSelect),2));
% [zValsRawAct, pValsRaw, actClust] = timePermCluster(ieegMeanSyllableRoi,ieegMeanSyllableShuffleRoi,1000,1,1.645);
 [zValsRawAct, pValsRaw, actClust] = timePermClusterAfterPerm(ieegMeanRoi,ieegShuffleRoi,1,1.645);
pmask = zeros(1,size(ieegMeanRoi,2));
% [~, ~, actClust]=timePermClusterZOptim(acctimeall,acctimeshuffle,10000,1,norminv(1-pVal2Cutoff));

for iClust=1:length(actClust.Size)
    if actClust.Size{iClust}>actClust.perm95
        pmask(actClust.Start{iClust}: ...
            actClust.Start{iClust}+(actClust.Size{iClust}-1)) ...
            = 1;
    end
end
visTimePlot(timeRange,ieegMeanRoiSubj,colval = colGroup(3,:), labels = xlab)
hold on;
visTimePlot(timeRange,squeeze((mean(ieegShuffMod(:,roiIds,:),1))),colval = [0.5 0.5 0.5], labels = xlab)
scatter(timeGammaMod(logical(pmask)),ylineval.*ones(1,sum(pmask)),'filled',MarkerEdgeColor=colGroup(3,:),MarkerFaceColor=colGroup(3,:));
ieegRoiEpochScaled = ieegRoi.*pmask; [maxAcc,maxId] = max(ieegRoiEpochScaled');
timeSyllableMaxGroup{3} = timeGammaMod(maxId);
time_centroid = [];
for iElec = 1:size(ieegRoiEpochScaled,1)
    time_centroid = [time_centroid calculateTimeCentroid(timeGammaMod,ieegRoiEpochScaled(iElec,:))];
   
end
time_centroid = time_centroid(time_centroid>timeStart);
timeSyllableCentroidGroup{3} = time_centroid;

ylim(yrange)

%%
    hemisphere = 'lh'

    bnLabels = extractBNlabels('triangular',hemisphere)
    roiIds = extractRoiIds(Subject, channelNamePooledSelect,bnLabels);
ieegMeanRoiTrial = squeeze(mean(ieegMeanMod(trials2select,roiIds,:),1));
[ieegMeanRoi,ieegMeanRoiSubj] = extractMeanControlSubject(ieegMeanRoiTrial,channelNamePooledSelect(roiIds));
 [~,ieegMeanShuffSubj] = extractMeanControlSubject(squeeze(mean(ieegShuffMod(:,roiIds,:),1)),channelNamePooledSelect(roiIds));
ieegShuffleRoi = [];
for iTer = 1:1000
    trials2select = randperm(size(ieegMeanMod,1),groupNum);
   % ieegShuffleRoi(iTer, :) = squeeze(mean(mean(ieegShuffMod(trials2select,roiIds,:),2),1))';
   ieegShuffRoiTrial = squeeze(mean(ieegShuffMod(trials2select,roiIds,:),1));
    ieegShuffleRoi(iTer, :) = extractMeanControlSubject(ieegShuffRoiTrial,channelNamePooledSelect(roiIds));
end
% ieegMeanSyllableRoi = [];
% ieegMeanSyllableShuffleRoi = [];
% 
% for iter = 1:500
%     [~,ieegMeanSyllableRoi(iter,:,:)] = extractMeanControlSubject( squeeze(ieegMeanMod(iter,roiIds,:)),channelNamePooledSelect(roiIds));
%     [~,ieegMeanSyllableShuffleRoi(iter,:,:)] = extractMeanControlSubject( squeeze(ieegShuffMod(iter,roiIds,:)),channelNamePooledSelect(roiIds));
% end
% ieegMeanSyllableRoi = squeeze(mean(ieegMeanSyllableRoi,1));
% ieegMeanSyllableShuffleRoi = squeeze(mean(ieegMeanSyllableShuffleRoi,1));
%  [zValsRawAct, pValsRaw, actClust] = timePermCluster(ieegMeanSyllableRoi,ieegMeanSyllableShuffleRoi,2000,1,1.645);
[zValsRawAct, pValsRaw, actClust] = timePermClusterAfterPerm(ieegMeanRoi,ieegShuffleRoi,1,1.645);
pmask = zeros(1,size(ieegMeanRoi,2));

for iClust=1:length(actClust.Size)
    if actClust.Size{iClust}>actClust.perm95
        pmask(actClust.Start{iClust}: ...
            actClust.Start{iClust}+(actClust.Size{iClust}-1)) ...
            = 1;
    end
end
figure;
visTimePlot(timeRange,ieegMeanRoiSubj,colval = colors(11,:), labels = xlab )
hold on;
visTimePlot(timeRange,squeeze((mean(ieegShuffMod(:,roiIds,:),1))),colval = [0.5 0.5 0.5], labels = xlab)
scatter(timeGammaMod(logical(pmask)),ylineval.*ones(1,sum(pmask)),'filled',MarkerEdgeColor=colors(11,:),MarkerFaceColor=colors(11,:));

title(['Pars-triangularis ' num2str(sum(roiIds))])

%%
    hemisphere = 'lh'
    
    bnLabels = extractBNlabels('heschl',hemisphere)
    hesch_roiIds = extractRoiIds(Subject, channelNamePooledSelect,bnLabels);
    heschlChannel = channelNamePooledSelect(hesch_roiIds);
    heshcSubj = extractSubjectPerChannel(heschlChannel)
    
    bnLabels = extractBNlabels('heschl',hemisphere)
    stg_roiIds = extractRoiIds(Subject, channelNamePooledSelect,bnLabels);
    stgChannel = channelNamePooledSelect(stg_roiIds);
    stg_heshSubj_ids = extractChannelPerSubject(stgChannel,heshcSubj);
    stg_heshSubj_chan = stgChannel(stg_heshSubj_ids);

    roiIds = ismember(channelNamePooledSelect,stg_heshSubj_chan);



% channelInfoAll = extractChannelLocation(Subject,channelNamePooledSelect);
% xyzChannel = reshape(extractfield(channelInfoAll,'xyz'),3,length(channelNamePooledSelect));
% preAcIds = xyzChannel(2,:)>0;
%roiIds = roiIds & preAcIds;

ieegRoi = squeeze(mean(ieegMeanMod(:,roiIds,:),1));

trials2select = randperm(size(ieegMeanMod,1),groupNum);
ieegMeanRoi = squeeze(mean(mean(ieegMeanMod(:,roiIds,:),2),1))';
ieegShuffleRoi = [];
for iTer = 1:1000
    trials2select = randperm(size(ieegMeanMod,1),groupNum);
    ieegShuffleRoi(iTer, :) = squeeze(mean(mean(ieegShuffMod(trials2select,roiIds,:),2),1))';
end
% [zValsRawAct, pValsRaw, actClust] = timePermCluster(ieegRoi,ieegShuffleRoi,1000,1,1.645);
 [zValsRawAct, pValsRaw, actClust] = timePermClusterAfterPerm(ieegMeanRoi,ieegShuffleRoi,1,1.645);
pmask = zeros(1,size(ieegMeanRoi,2));
% [~, ~, actClust]=timePermClusterZOptim(acctimeall,acctimeshuffle,10000,1,norminv(1-pVal2Cutoff));

for iClust=1:length(actClust.Size)
    if actClust.Size{iClust}>actClust.perm95
        pmask(actClust.Start{iClust}: ...
            actClust.Start{iClust}+(actClust.Size{iClust}-1)) ...
            = 1;
    end
end
figure;
visTimePlot(timeRange,squeeze((mean(ieegMeanMod(:,roiIds,:),1))),colval = colors(11,:), labels = xlab )
hold on;
visTimePlot(timeRange,squeeze((mean(ieegShuffMod(:,roiIds,:),1))),colval = [0.5 0.5 0.5], labels = xlab)
scatter(timeGammaMod(logical(pmask)),ylineval.*ones(1,sum(pmask)),'filled',MarkerEdgeColor=colors(11,:),MarkerFaceColor=colors(11,:));
groupingIds(roiIds) = 7;
title(' Heschl')
ylim([-0.01 0.3])
%% 

ieegMeanModPeak = max(squeeze(mean(ieegMeanMod,1))');
plot_activation_density_v3(channelNamePooledSelect, ieegMeanModPeak', cLim = [0 0.2],...
    gaussFwhm=20, diskThresh = 100, stype='pial', ptype='inflated', colbarTitle =' HG power',...
    transparentPoint=0.1, colorscale = 'hot');

%% ROI histogram plot
labels = {'IFG','cMFG','rMFG','SMC','IPC','Insula','aSTG','pSTG','STS','PAC'};
roiColors = [colors(4,:); colors(3,:); colors(10,:); colors(5,:); colors(2,:);...
    colors(6,:); colors(9,:); colors(1,:); colors(8,:); colors(7,:);];
% medX = cellfun(@median, timeCentroid);
% 
% % Sort medians and get the sort order
% [~, sortIdx] = sort(medX);
    
timeMaxRoi2plot = timeSyllableCentroid; 
labels2plot = labels;

plotRoiHistograms(timeMaxRoi2plot,labels2plot,roiColors)


%% Group histogram plot

% timeGroup{1} = [timeSyllableCentroid{1:3}];
% timeGroup{2} = [timeSyllableCentroid{4:5}];
% timeGroup{3} = [timeSyllableCentroid{6:10}];

labels2plot = {'planning','articulation','monitoring'};
 plotRoiHistograms(timeSyllableCentroidGroup,labels2plot,lines(3))

 [h,p,ci,stats] = ttest2( timeSyllableCentroidGroup{2},timeSyllableCentroidGroup{1})
 [h,p,ci,stats] = ttest2( timeSyllableCentroidGroup{3},timeSyllableCentroidGroup{2})
 [h,p,ci,stats] = ttest2( timeSyllableCentroidGroup{3},timeSyllableCentroidGroup{1})

%% 

% Example cell array with three groups
data = timeSyllableCentroidGroup;

% Convert cell array data to a format suitable for ANOVA
group = [];
values = [];

for i = 1:numel(data)
    group = [group; i*ones(size(data{i},2), 1)]; % Assign group numbers
    values = [values; data{i}']; % Concatenate all values into a single array
end
group3 = group;
group3(group<=3) = 1;
group3(group>3&group<=5) = 2;
group3(group>5) = 3;

% Perform one-way ANOVA
[p, tbl, stats] = anova1(values, group, 'off');  % 'off' turns off the ANOVA table display

% Display results
disp('ANOVA Table:');
disp(tbl);
fprintf('p-value = %.4f\n', p);

% If you want to investigate significant differences further with post-hoc tests
if p < 0.05
    figure;
    [c,m,h,gnames] = multcompare(stats);
    title('Multiple Comparisons of Means');
end

%%
accChanTime = [];
accChanTimeChance = [];
for iChan = 1:128
    for iPos = 1:3
        accTime = zeros(1,281);
        accTimeChance = zeros(1,281);
        for iter =1:10
            accTime = accTime + decodeTimeStruct1D{iChan,iPos,iter}.accTime;
            accTimeChance = accTimeChance + decodeTimeStruct1Dshuffle{iChan,iPos,iter}.accTime;
        end
        accTime = accTime./10;
        accTimeChance = accTimeChance./10;
        accChanTime(iChan,iPos,:) = accTime;
        accChanTimeChance(iChan,iPos,:) = accTimeChance;
    end
end

figure;
for iChan = 1:128
    subplot(size(chanMap,1),size(chanMap,2),find(ismember(chanMap',iChan)));
     plot(squeeze(accChanTimeChance(iChan,:,:))')
     ylim([0 0.3])
end

figure;
for iChan = 1:128
    subplot(size(chanMap,1),size(chanMap,2),find(ismember(chanMap',iChan)));
     plot(squeeze(accChanTime(iChan,:,:)-accChanTimeChance(iChan,:,:))')
     ylim([-0.1 0.2])
end
%% Subject specific analysis - power traces
subjNames = {'D35','D39','D41','D49','D59','D64','D68','D79','D81','D88','D93','D96'};
%subjNames = extractSubjectPerChannel(channelNameEpoch);    
for iSubj = 1:length(subjNames)
    iSubj
    channelIds = extractChannelPerSubject(channelNameEpoch,subjNames{iSubj});
    channelNameEpoch(channelIds)

    ieegEpochSubj = ieegEpoch(channelIds,:);
    ieegBaseSubj = ieegBase(channelIds,:);
    figure;
    bnLabels = [extractBNlabels('ifg',hemisphere) extractBNlabels('rmfg',hemisphere) extractBNlabels('cmfg',hemisphere) ];
    roiIds = extractRoiIds(Subject, channelNameEpoch(channelIds),bnLabels);    
    if(sum(roiIds)>1)
    ieegRoiEpoch = ieegEpochSubj(roiIds,:);
    ieegRoiBase = ieegBaseSubj(roiIds,:);
    
    pmask = zeros(1,size(ieegRoiBase,2));
    [zValsRawAct, pValsRaw, actClust] = timePermCluster(ieegRoiEpoch,ieegRoiBase,1000,1,1.645);
    
    for iClust=1:length(actClust.Size)
        if actClust.Size{iClust}>actClust.perm95
            pmask(actClust.Start{iClust}: ...
                actClust.Start{iClust}+(actClust.Size{iClust}-1)) ...
                = 1;
        end
    end    
    ieegRoiEpochScaled = ieegRoiEpoch.*pmask;

    timeCentroidRoi = [];
    for iElec = 1:size(ieegRoiEpochScaled,1)
        timeCentroidRoi = [timeCentroidRoi calculateTimeCentroid(timeGamma,ieegRoiEpochScaled(iElec,:))];   
        
    end
    timeCentroidRoi = timeCentroidRoi(timeCentroidRoi>-0.9);
    timeCentroid{iSubj,1} = timeCentroidRoi;
    
    visTimePlot(timeRange,ieegRoiEpoch,colval = [0.05 0.03 0.53], labels = xlab )
    hold on;
    yline(0, '-','','LineWidth',1, 'Color','k');

    scatter(timeGamma(logical(pmask)),-0.1.*ones(1,sum(pmask)),'filled',...
    MarkerEdgeColor=[0.05 0.03 0.53],MarkerFaceColor=[0.05 0.03 0.53]);
    end
    bnLabels = [extractBNlabels('smc',hemisphere) extractBNlabels('ipc',hemisphere)];
    roiIds = extractRoiIds(Subject, channelNameEpoch(channelIds),bnLabels);    
    if(sum(roiIds)>1)
    ieegRoiEpoch = ieegEpochSubj(roiIds,:);
    ieegRoiBase = ieegBaseSubj(roiIds,:);
    
    pmask = zeros(1,size(ieegRoiBase,2));
    [zValsRawAct, pValsRaw, actClust] = timePermCluster(ieegRoiEpoch,ieegRoiBase,1000,1,1.645);
    
    for iClust=1:length(actClust.Size)
        if actClust.Size{iClust}>actClust.perm95
            pmask(actClust.Start{iClust}: ...
                actClust.Start{iClust}+(actClust.Size{iClust}-1)) ...
                = 1;
        end
    end    
    ieegRoiEpochScaled = ieegRoiEpoch.*pmask;
    timeCentroidRoi = [];
    for iElec = 1:size(ieegRoiEpochScaled,1)
        timeCentroidRoi = [timeCentroidRoi calculateTimeCentroid(timeGamma,ieegRoiEpochScaled(iElec,:))];   
        
    end
    timeCentroidRoi = timeCentroidRoi(timeCentroidRoi>-0.9);
    timeCentroid{iSubj,2} = timeCentroidRoi;
    visTimePlot(timeRange,ieegRoiEpoch,colval = [0.8 0.28 0.47], labels = xlab )
    hold on;
    yline(0, '-','','LineWidth',1, 'Color','k');

    scatter(timeGamma(logical(pmask)),-0.2.*ones(1,sum(pmask)),'filled',...
    MarkerEdgeColor=[0.8 0.28 0.47],MarkerFaceColor=[0.8 0.28 0.47]);
    end

    bnLabels = [extractBNlabels('astg',hemisphere) extractBNlabels('pstg',hemisphere) extractBNlabels('insula',hemisphere) extractBNlabels('heschl',hemisphere) extractBNlabels('sts',hemisphere)];
    roiIds = extractRoiIds(Subject, channelNameEpoch(channelIds),bnLabels);    
    if(sum(roiIds)>1)
    ieegRoiEpoch = ieegEpochSubj(roiIds,:);
    ieegRoiBase = ieegBaseSubj(roiIds,:);
    
    pmask = zeros(1,size(ieegRoiBase,2));
    [zValsRawAct, pValsRaw, actClust] = timePermCluster(ieegRoiEpoch,ieegRoiBase,1000,1,1.645);
    
    for iClust=1:length(actClust.Size)
        if actClust.Size{iClust}>actClust.perm95
            pmask(actClust.Start{iClust}: ...
                actClust.Start{iClust}+(actClust.Size{iClust}-1)) ...
                = 1;
        end
    end    
    ieegRoiEpochScaled = ieegRoiEpoch.*pmask;
    timeCentroidRoi = [];
    for iElec = 1:size(ieegRoiEpochScaled,1)
        timeCentroidRoi = [timeCentroidRoi calculateTimeCentroid(timeGamma,ieegRoiEpochScaled(iElec,:))];   
        
    end
    timeCentroidRoi = timeCentroidRoi(timeCentroidRoi>-0.9);
    timeCentroid{iSubj,3} = timeCentroidRoi;
    visTimePlot(timeRange,ieegRoiEpoch,colval = [0.13 0.57 0.55], labels = xlab )
    hold on;
    yline(0, '-','','LineWidth',1, 'Color','k');

    scatter(timeGamma(logical(pmask)),-0.3.*ones(1,sum(pmask)),'filled',...
    MarkerEdgeColor=[0.13 0.57 0.55],MarkerFaceColor=[0.13 0.57 0.55]);
    end
    title(subjNames{iSubj})
end
timeCentMedian = [];
for iSubj = 1:length(subjNames)
%     figure;
%     labels2plot = {'planning','articulation','monitoring'};
%     plotRoiHistograms( timeCentroid(iSubj,:),labels2plot,[0.05 0.03 0.53; 0.8 0.28 0.47; 0.13 0.57 0.55])
%     sgtitle(subjNames{iSubj})
    timeCentMedian(iSubj,1) = median(timeCentroid{iSubj,1});
    if(timeCentMedian(iSubj,1)>1)
        timeCentMedian(iSubj,1) = nan;
    end
    timeCentMedian(iSubj,2) = median(timeCentroid{iSubj,2});
    timeCentMedian(iSubj,3) = median(timeCentroid{iSubj,3});

end

figure;
scatter(timeCentMedian(:,1)', 1:12);
hold on;
scatter(timeCentMedian(:,2)', 1:12);
scatter(timeCentMedian(:,3)', 1:12);

%% Subject specific analysis - syllable

subjNames = {'D35','D39','D41','D49','D59','D64','D68','D79','D81','D88','D93','D96'};
%subjNames = extractSubjectPerChannel(channelNamePooled);    
for iSubj = 1:length(subjNames)
    iSubj
    channelIds = extractChannelPerSubject(channelNamePooled,subjNames{iSubj});
    channelNamePooled(channelIds)

    ieegMeanSubj = ieegMeanMod(:,channelIds,:);
    ieegShuffSubj = ieegShuffMod(:,channelIds,:);


    bnLabels = [extractBNlabels('ifg',hemisphere) extractBNlabels('rmfg',hemisphere) extractBNlabels('cmfg',hemisphere) ];
   
    roiIds = extractRoiIds(Subject, channelNamePooled(channelIds),bnLabels);

    if(sum(roiIds)>1)
    
    trials2select = randperm(size(ieegMeanMod,1),groupNum);
    ieegMeanRoi = squeeze(mean(mean(ieegMeanSubj(trials2select,roiIds,:),2),1))';
    
    ieegShuffleRoi = [];
    for iTer = 1:1000
        trials2select = randperm(size(ieegMeanMod,1),groupNum);
        ieegShuffleRoi(iTer, :) = squeeze(mean(mean(ieegShuffSubj(trials2select,roiIds,:),2),1))';       
    end
    
     [zValsRawAct, pValsRaw, actClust] = timePermClusterAfterPerm(ieegMeanRoi,ieegShuffleRoi,1,1.645);
    pmask = zeros(1,size(ieegMeanRoi,2));
    % [~, ~, actClust]=timePermClusterZOptim(acctimeall,acctimeshuffle,10000,1,norminv(1-pVal2Cutoff));
    
    for iClust=1:length(actClust.Size)
        if actClust.Size{iClust}>actClust.perm95
            pmask(actClust.Start{iClust}: ...
                actClust.Start{iClust}+(actClust.Size{iClust}-1)) ...
                = 1;
        end
    end
    ieegRoi = squeeze(mean(ieegMeanMod(:,roiIds,:),1));
    ieegRoiEpochScaled = ieegRoi.*pmask;
    timeCentroidRoi = [];
    for iElec = 1:size(ieegRoiEpochScaled,1)
        timeCentroidRoi = [timeCentroidRoi calculateTimeCentroid(timeGammaMod,ieegRoiEpochScaled(iElec,:))];   
        
    end
    timeCentroidRoi = timeCentroidRoi(timeCentroidRoi>-0.9);
    timeCentroid{iSubj,1} = timeCentroidRoi;
    
    figure;
    visTimePlot(timeRange,squeeze((mean(ieegMeanSubj(:,roiIds,:),1))),colval = colors(1,:), labels = xlab )
    hold on;
    visTimePlot(timeRange,squeeze((mean(ieegShuffSubj(:,roiIds,:),1))),colval = [0.5 0.5 0.5], labels = xlab)
    scatter(timeGammaMod(logical(pmask)),-0.005.*ones(1,sum(pmask)),'filled',MarkerEdgeColor=colors(1,:),MarkerFaceColor=colors(1,:));
    end

    bnLabels = [extractBNlabels('smc',hemisphere) extractBNlabels('ipc',hemisphere)];
    roiIds = extractRoiIds(Subject, channelNamePooled(channelIds),bnLabels);
    if(sum(roiIds)>1)
    
    
    trials2select = randperm(size(ieegMeanMod,1),groupNum);
    ieegMeanRoi = squeeze(mean(mean(ieegMeanSubj(trials2select,roiIds,:),2),1))';
    
    ieegShuffleRoi = [];
    for iTer = 1:1000
        trials2select = randperm(size(ieegMeanMod,1),groupNum);
        ieegShuffleRoi(iTer, :) = squeeze(mean(mean(ieegShuffSubj(trials2select,roiIds,:),2),1))';       
    end
    
     [zValsRawAct, pValsRaw, actClust] = timePermClusterAfterPerm(ieegMeanRoi,ieegShuffleRoi,1,1.645);
    pmask = zeros(1,size(ieegMeanRoi,2));
    % [~, ~, actClust]=timePermClusterZOptim(acctimeall,acctimeshuffle,10000,1,norminv(1-pVal2Cutoff));
    
    for iClust=1:length(actClust.Size)
        if actClust.Size{iClust}>actClust.perm95
            pmask(actClust.Start{iClust}: ...
                actClust.Start{iClust}+(actClust.Size{iClust}-1)) ...
                = 1;
        end
    end
    ieegRoi = squeeze(mean(ieegMeanMod(:,roiIds,:),1));
    ieegRoiEpochScaled = ieegRoi.*pmask;
    timeCentroidRoi = [];
    for iElec = 1:size(ieegRoiEpochScaled,1)
        timeCentroidRoi = [timeCentroidRoi calculateTimeCentroid(timeGammaMod,ieegRoiEpochScaled(iElec,:))];   
        
    end
    timeCentroidRoi = timeCentroidRoi(timeCentroidRoi>-0.9);
    timeCentroid{iSubj,2} = timeCentroidRoi;
    
    visTimePlot(timeRange,squeeze((mean(ieegMeanSubj(:,roiIds,:),1))),colval = colors(2,:), labels = xlab )
    hold on;
    visTimePlot(timeRange,squeeze((mean(ieegShuffSubj(:,roiIds,:),1))),colval = [0.5 0.5 0.5], labels = xlab)
    scatter(timeGammaMod(logical(pmask)),-0.01.*ones(1,sum(pmask)),'filled',MarkerEdgeColor=colors(2,:),MarkerFaceColor=colors(2,:));
    end


    bnLabels = [extractBNlabels('astg',hemisphere) extractBNlabels('pstg',hemisphere) extractBNlabels('insula',hemisphere) extractBNlabels('heschl',hemisphere) extractBNlabels('sts',hemisphere)];
    roiIds = extractRoiIds(Subject, channelNamePooled(channelIds),bnLabels);
    if(sum(roiIds)>1)
    
    
    trials2select = randperm(size(ieegMeanMod,1),groupNum);
    ieegMeanRoi = squeeze(mean(mean(ieegMeanSubj(trials2select,roiIds,:),2),1))';
    
    ieegShuffleRoi = [];
    for iTer = 1:1000
        trials2select = randperm(size(ieegMeanMod,1),groupNum);
        ieegShuffleRoi(iTer, :) = squeeze(mean(mean(ieegShuffSubj(trials2select,roiIds,:),2),1))';       
    end
    
     [zValsRawAct, pValsRaw, actClust] = timePermClusterAfterPerm(ieegMeanRoi,ieegShuffleRoi,1,1.645);
    pmask = zeros(1,size(ieegMeanRoi,2));
    % [~, ~, actClust]=timePermClusterZOptim(acctimeall,acctimeshuffle,10000,1,norminv(1-pVal2Cutoff));
    
    for iClust=1:length(actClust.Size)
        if actClust.Size{iClust}>actClust.perm95
            pmask(actClust.Start{iClust}: ...
                actClust.Start{iClust}+(actClust.Size{iClust}-1)) ...
                = 1;
        end
    end
    ieegRoi = squeeze(mean(ieegMeanMod(:,roiIds,:),1));
    ieegRoiEpochScaled = ieegRoi.*pmask;
    timeCentroidRoi = [];
    for iElec = 1:size(ieegRoiEpochScaled,1)
        timeCentroidRoi = [timeCentroidRoi calculateTimeCentroid(timeGammaMod,ieegRoiEpochScaled(iElec,:))];   
        
    end
    timeCentroidRoi = timeCentroidRoi(timeCentroidRoi>-0.9);
    timeCentroid{iSubj,3} = timeCentroidRoi;
    
    visTimePlot(timeRange,squeeze((mean(ieegMeanSubj(:,roiIds,:),1))),colval = colors(3,:), labels = xlab )
    hold on;
    visTimePlot(timeRange,squeeze((mean(ieegShuffSubj(:,roiIds,:),1))),colval = [0.5 0.5 0.5], labels = xlab)
    scatter(timeGammaMod(logical(pmask)),-0.015.*ones(1,sum(pmask)),'filled',MarkerEdgeColor=colors(3,:),MarkerFaceColor=colors(3,:));
    end
    title(subjNames{iSubj})

end

