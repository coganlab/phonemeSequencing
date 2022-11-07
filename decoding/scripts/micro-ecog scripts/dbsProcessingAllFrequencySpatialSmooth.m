addpath(genpath('C:\Users\sd355\Box\ECoG_Recon\matlab_code\'));
global DUKEDIR
DUKEDIR = 'C:\Users\sd355\Box\CoganLab\Data\Micro\Processed Data';
dLabels = dir(DUKEDIR);
dLabels = dLabels(3:end);

tw = [-2 2]; % time window


baseTW = [-0.5 0]; % preonset time window
activeTW = [-1 1]; % postonset time window
gammaF = [70 150]; % frequency in Hz
fsDown = 200;
load('channelMap.mat');
chanMap = (chanMap');
selectedChannels = sort(chanMap(~isnan(chanMap)))';
%% data Loading
subjectId = 'S14';

Experiment = loadExperiment(subjectId);
fsD = Experiment.processing.ieeg.sample_rate;
Trials = dbTrials(subjectId,Experiment.recording.recording_day,'Speech_OvertMimeMove');

trialFiles = strcat('\',Experiment.recording.recording_day,'\mat\trialInfo.mat');
[ieegResponse,~,trigOnset]=trialIEEGUpdate(Trials,selectedChannels,'phon1Onset','ieeg',tw.*1000);
[ieegAuditory,~,trigAuditory]=trialIEEGUpdate(Trials,selectedChannels,'Auditory','ieeg',tw.*1000);

ieegAuditory = permute(ieegAuditory,[2,1,3]);
ieegResponse = permute(ieegResponse,[2,1,3]);
 
respId = find(~isnan(trigOnset));
ieegResponse = ieegResponse(:,respId,:);
ieegAuditory = ieegAuditory(:,respId,:); 
Trials = Trials(respId);
  
trigOns = trigOnset(respId)./fsD;
trigAud = trigAuditory(respId)./fsD;
respTime  = trigOns-trigAud;
[respTimeSort,sortId] = sort(respTime);
%% Ieeg Class definition

ieegBase = ieegStructMicro(ieegAuditory, fsD, tw, [1 fsD/2], 'Auditory', chanMap);
ieegResponse = ieegStructMicro(ieegResponse, fsD, tw, [1 fsD/2], 'Response', chanMap);

%% Phoneme trial parser
phonemeTrial = phonemeSequenceTrialParser(Trials);
%% Common average referencing
load('impedance.mat')
higImpId = find(log10(impedance1)>6);
ieegBaseCar = extractCar(ieegBase, higImpId);
clear ieegBase
ieegResponseCar = extractCar(ieegResponse, higImpId);
clear ieegResponse
%% Extract High Gamma

ieegHiGammaNormResp = extractHiGammaNorm(ieegResponseCar,ieegBaseCar,fsDown,activeTW,baseTW);

%% Phoneme results

decoderClass = phonemeDecoderClass(5,80);
decodeResultStruct = baseClassify(decoderClass,ieegHiGammaNormResp,phonemeTrial,'phoneme',[-0.5 0.5]);

%% Phoneme maps
decodeResultStruct = indChanClassify(decoderClass,ieegHiGammaNormResp,phonemeTrial,'class',[-0.5 0.5]);
chanLabAuc = [];
for iChan = 1:length(decodeResultStruct)
    chanLabAuc(iChan,:) = decodeResultStruct{iChan}.aucAll(1,:); 
end
labels = {'low','high','labials','dorsals'};
%labels = {'a','ae','i','u','b','p','v','g','k'};

%chanLabDiff = abs(chanLabAuc - mean(chanLabAuc,2));
chanLabDiff = (chanLabAuc - 0.5);
figure;

for iPhon = 1:4
     subplot(2,2,iPhon)
%figure
    chanVal = chanView(squeeze(chanLabAuc(:,iPhon)),(chanMap),selectedChannels,[],labels{iPhon},[0.6 1],[],0.25);
    
    colormap(flipud(bone(4096)));
    prc4cutoff = 90;
    cutoff = prctile(chanVal(:),prc4cutoff)
    hold on;
    
    
    chanValMark = true(size(chanVal));
    chanValMark(chanVal<cutoff) = false;
    chanValTemp = chanVal;
    chanValTemp(chanVal<cutoff) = 0;
    measurements = regionprops(chanValMark, chanValTemp, 'all');
    
    
   
    set(gca,'FontSize',15);
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    axis equal
    axis tight
    if(isempty(measurements)||round(cutoff,2)<0.75)
        continue;
    else
     [intensePhoneme(iPhon),maxId] = max([measurements(:).MaxIntensity]);
    hold on;
    scatter(measurements(maxId).ConvexHull(:,1)',measurements(maxId).ConvexHull(:,2)',20,'r','filled')
    x = sqrt(measurements(maxId).ConvexArea);
    phonArea = (1.33*x)^2
    end
    centerOfMass = measurements(maxId).WeightedCentroid;
    scatter(centerOfMass(1),centerOfMass(2),50,'x','LineWidth',2)
    
end

%%
elecComb = nchoosek(1:length(selectedChannels),2);
eucDist = [];
pitch = 1.33;
sigCov = [];
for eD = 1 : size(elecComb,1)
    eD
    [x1,y1] = find(ismember(chanMap,selectedChannels(elecComb(eD,1))));
    [x2,y2] = find(ismember(chanMap,selectedChannels(elecComb(eD,2))));
    eucDist(eD) = pitch.*sqrt((x2-x1)^2+(y2-y1)^2);
    for iPhoneme = 1:4
        gamma1 = squeeze(chanLabAuc((elecComb(eD,1)),iPhoneme));
        gamma2 = squeeze(chanLabAuc((elecComb(eD,2)),iPhoneme));
        covtemp = cov(gamma1,gamma2);
        sigCov(eD,iPhoneme) = sqrt(-2*gamma1*gamma2+(gamma1^2+gamma2^2));
        %sigCov(eD,iPhoneme) = sqrt(gamma1*gamma2);
        
    end
end

 for iPhoneme = 1:4
    eucDistUnique = unique(eucDist);
    for eDU = 1:length(eucDistUnique)
        eIDs = find(eucDist == eucDistUnique(eDU));
        sigCovMean(eDU,iPhoneme) = median(sigCov(eIDs,iPhoneme));
        sigCovStd(eDU,iPhoneme) = std(sigCov(eIDs,iPhoneme))./sqrt(length(eIDs));   
    end
 end
% figure;
% errorbar(eucDistUnique, sigCovMean(:,2)',sigCovStd(:,2)');
% hold on;
% errorbar(eucDistUnique, sigCovMean(:,7)',sigCovStd(:,7)');
% errorbar(eucDistUnique, sigCovMean(:,11)',sigCovStd(:,11)');
% errorbar(eucDistUnique, sigCovMean(:,12)',sigCovStd(:,12)');

figure;
for iPhoneme = 1:4
errorbar(eucDistUnique, sigCovMean(:,iPhoneme)',sigCovStd(:,iPhoneme)');
hold on;
end

axis square
xlabel('Distance (mm)');
ylabel('Semivariance');
set(gca,'FontSize',20);
xlim([0 10])
legend('low','high','labial','dorsal');

%% Extract Spatial Smooth
accSmooth = zeros(10,8);
for iWindow = 1:8
    iWindow
    if(iWindow == 1)
        ieegHiGammaNormResp = extractHiGammaNorm(ieegResponseCar,ieegBaseCar,200,[-0.5 0.5],[-0.5 0]);
        
    else
        ieegBaseCarSmooth = spatialSmoothMicro(ieegBaseCar,[iWindow iWindow]);
        ieegResponseCarSmooth = spatialSmoothMicro(ieegResponseCar,[iWindow iWindow]);
        ieegHiGammaNormResp = extractHiGammaNorm(ieegResponseCarSmooth,ieegBaseCarSmooth,200,[-0.5 0.5],[-0.5 0]);   
        
    end
    for iTer = 1:10
        decoderClass = phonemeDecoderClass(20,80);
        decodeResultStruct = baseDecoder(decoderClass,ieegHiGammaNormResp,phonemeTrial,'phoneme',[-0.25 0.25]);
        accSmooth(iTer,iWindow) = decodeResultStruct.accPhonemeUnBias(1);
    end
end



elecSpaceSampStr = [];
accSpaceSampAll = [];
phonErrorSpaceSampAll = [];
for iWindow = 1:8  
    elecSpaceSampStrTemp = [num2str(iWindow) 'x' num2str(iWindow)];
    elecSpaceSampStr = [elecSpaceSampStr; repmat({elecSpaceSampStrTemp},10,1)];
    accSpaceSampAll = [accSpaceSampAll accSmooth(:,iWindow)'];
    accSpaceMean(iWindow) = median(accSmooth(:,iWindow)');
   
end 
save('PhonemeDecodeSpaceSmoothFirstPhoneme.mat','accSpaceSampAll','phonErrorSpaceSampAll','elecSpaceSampStr');

 h = figure('Units', 'pixels', ...
    'Position', [100 100 750 500]);
set(h, 'DefaultTextFontSize', 15); h = boxplot(accSpaceSampAll.*100,elecSpaceSampStr','symbol','','Colors','k');
set(h,{'linew'},{2})
ylim([0 0.7.*100]);
hold on;
plot(1:8,accSpaceMean.*100,'LineWidth',2);
yline(0.1111.*100, '--','chance','LineWidth',2, 'LabelHorizontalAlignment','left');
yline(max(accSpaceMean.*100).*0.95, ':','95% threshold','LineWidth',2, 'Color','r', 'LabelHorizontalAlignment','left');
xlabel('Spatial-average grid size');
ylabel('P1 accuracy (%)');
set(gca,'FontSize',15);
axis square
set(gca,'FontSize',20);


axis square

set( gca                       , ...
    'FontName'   , 'Arial' );


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
%% Phoneme results
numFold = 20;
varExplained = 80;
decoderClass = phonemeDecoderClass(numFold,varExplained);
decodeResultStruct = baseDecoder(decoderClass,ieegHiGammaNormResp,phonemeTrial,'phoneme',[-0.25 0.25]);

%% Tsne score
decompStructPhoneme = tsneDecompose(decoderClass,ieegHiGammaNormResp,phonemeTrial,'phoneme', [-0.25 0.25],[512 120 60 30],50);
decompStructClass = tsneDecompose(decoderClass,ieegHiGammaNormResp,phonemeTrial,'class', [-0.25 0.25],[512 120 60 30],50);

figure;
subplot(2,1,1);
h = boxplot(decompStructClass.silScore', 'symbol','');
set(h,{'linew'},{2})
axis square
set(gca,'XTickLabels', {'100','50','25','12.5'});
set(gca,'FontSize',10);
ylabel('Silhouette score');
title('Articulator grouping');
ylim([0 max(decompStructClass.silScore(:))]);
subplot(2,1,2);
h = boxplot(decompStructPhoneme.silScore', 'symbol','');
set(h,{'linew'},{2})
axis square
set(gca,'XTickLabels', {'100','50','25','12.5'});
set(gca,'FontSize',10);
xlabel('% Sampling');
title('Phoneme grouping');
ylim([0 max(decompStructPhoneme.silScore(:))]);

figure;
subplot(2,1,1);
h = boxplot(decompStructClass.silRatio', 'symbol','');
set(h,{'linew'},{2})
axis square
set(gca,'XTickLabels', {'100','50','25','12.5'});
set(gca,'FontSize',10);
ylabel('Silhouette ratio');
title('Phoneme grouping');
ylim([0 max(decompStructClass.silRatio(:))]);
subplot(2,1,2);
h = boxplot(decompStructPhoneme.silRatio', 'symbol','');
set(h,{'linew'},{2})
axis square
set(gca,'XTickLabels', {'100','50','25','12.5'});
set(gca,'FontSize',10);
xlabel('% Sampling');
title('Articulator grouping');
ylim([0 max(decompStructPhoneme.silRatio(:))]);
