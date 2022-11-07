addpath(genpath('C:\Users\sd355\Box\CoganLab\matlab_code\'));
global DUKEDIR
DUKEDIR = 'C:\Users\sd355\Box\CoganLab\Data\Micro\Processed Data';

tw = [-2 2]; % time window
etw = [-1.5 1.5];
etwG = [-1 1];% epoch time window
%% data Loading
Experiment = loadExperiment('S14');
fsD = Experiment.processing.ieeg.sample_rate;
chanMap = Experiment.recording.channel_map;
chanMap = (chanMap');
selectedChannels = sort(chanMap(~isnan(chanMap)))';



fullPathToRawData = 'C:\Users\sd355\Box\CoganLab\uECoG_Upload\S14_IntraOp_Recording_12_03_19\';
fileId = 5:35;
fs = 20000;
[ieegall,~,trigger] = intanWrap(fullPathToRawData,selectedChannels,fileId,0,1);
timeAll = linspace(0,(size(ieegall,2)-1)./fs,size(ieegall,2));


Trials = dbTrials('S14',Experiment.recording.recording_day,'Speech_OvertMimeMove');
trialFiles = strcat('\',Experiment.recording.recording_day,'\mat\trialInfo.mat');
trigSpeech = [Trials.ResponseOnset]./30e3;

 

 phon1 = [Trials.phon1];


respId = find(~isnan(trigSpeech));
trigOns = trigSpeech(respId);
ieegSplit = splitIeeg(ieegall,round(trigOns.*fs),tw,fs);

phon1 = phon1(respId);
timeSplit = linspace(tw(1),tw(2),size(ieegSplit,3));

%% Common average referencing
% higImpId = log10(impedance1)>6;
% ieegCarImp = carFilterImpedance(ieegSplit,higImpId);

%% Spectrogram analysis
% [specRaw,pPerc,pProd]= getSpectrograms(ieegSplit,[],tw,etwG,[1 2000],[-0.75 -0.5],[-0.25 0.25],[1 1.5],[200 1500],fs,1);
% [~, pvalsMCleanPerc] = fdr(pPerc,0.01);
%     [~, pvalsMCleanProd] = fdr(pProd,0.05);
% figure;
% specChanMap(specRaw,(chanMap),1:128,find(pvalsMCleanPerc),etwG,[-0.75 -0.5],[0.5 1],[1 2000],[300 1750],[0.7 1.4],0);

%% Spike band filtering
ieegSpikeBand = [];


for iChan = 1:size(ieegSplit,1)
ieegSpikeBand(iChan,:,:) = eegfilt(squeeze(ieegSplit(iChan,:,:)),fs,300,2000,0,200);

end
% [~,goodtrials] = remove_bad_trials(ieegSplit,14);
% goodTrialsCommon = extractCommonTrials(goodtrials);
%% Spike detection and extraction

spikeData = [];

for iChan = 1:128
    iChan
    ieegSpikeBandTemp = squeeze(ieegSpikeBand(iChan,:,:));
    sizeIeeg = size(squeeze(ieegSpikeBandTemp));
    spikeRast = zeros(size(squeeze(ieegSpikeBandTemp)));
    spikeSplitChan = [];
    timeSpikeChan = [];
    spikeThresh = 4*rms((ieegSpikeBandTemp(:)));
    for iTrial = 1:size(ieegSpikeBandTemp,1)
        
        %spikeThresh = 3.75*rms(abs(ieegSpikeBandTemp(iTrial,:)));
       % spikeId = find(abs(ieegSpikeBandTemp(iTrial,:))>spikeThresh);  
         spikeId = spike_detection(abs(ieegSpikeBandTemp(iTrial,:)),spikeThresh);
         spikeId(diff(spikeId)<100)=[];
          spikeId(spikeId<100) = [];
          spikeId(spikeId+100>sizeIeeg(2)) = [];
         if(isempty(spikeId))
            continue
         else
             spikeSplitTrial = squeeze(splitIeeg(ieegSpikeBandTemp(iTrial,:),spikeId,[-0.0025 0.0025],fs));
             if(length(spikeId)==1)
                 spikeSplitChan = [spikeSplitChan; spikeSplitTrial']; 
                 timeSpikeChan = [timeSpikeChan timeSplit(spikeId)];
             else
                spikeSplitChan = [spikeSplitChan; spikeSplitTrial]; 
                timeSpikeChan = [timeSpikeChan timeSplit(spikeId)];
             end
%               figure; 
%                 plot(timeSplit,ieegSpikeBandTemp(iTrial,:));
%                 hold on;
%                 scatter(timeSplit(spikeId),spikeThresh.*ones(1,length(spikeId)),'filled');
         end
%        
        spikeRast(iTrial,spikeId)=1;
    end
%     figure;
%     timeSpike = linspace(-0.0025, 0.0025, size(spikeSplitChan,2));
%      plot(timeSpike,spikeSplitChan')
%     figure;
%     plotSpikeRaster(logical(spikeRast),'PlotType','scatter','XLim',[-1 1]);
    spikeData{iChan}.spikeRast = spikeRast;

    % PCA on spikes to see if they are positive are negative
    [~,score] = pca(spikeSplitChan);
    posSpike = score(:,1)>10;
    spikeData{iChan}.spikeSplitPos = spikeSplitChan(posSpike,:);
    spikeData{iChan}.spikeSplitNeg = spikeSplitChan(~posSpike,:);
end

%save('surfaceSpikes-Threshold3.75.mat','spikeSplitPos','spikeSplitNeg','spikeRastChan','-v7.3');
%% Channel raster plot

timeLab = linspace(-2,2,5);
spikeSmooth = [];
figure;
for iChan = 1:128
    spikeRast = spikeData{iChan}.spikeRast;
%     gaussFilt = gausswin(100);
%     
%     for iTrial = 1:size(spikeRast,1)
%         spikeSmooth(iChan, iTrial,:) = resample(filter(gaussFilt,1,spikeRast(iTrial,:)),200,fs);
%     end
%     figure;
%     imagesc(spikeSmooth);
    subplot(size(chanMap,1),size(chanMap,2),find(ismember(chanMap',iChan)));
   %  figure;
    markerFormat = struct();
    markerFormat.LineWidth =0.05;
    LineFormat = struct();
    LineFormat.LineWidth = 3;
     plotSpikeRaster(logical(spikeRast),'PlotType','scatter', 'MarkerFormat',markerFormat,  'LineFormat',LineFormat,'XLimForCell',[-2 2]);
%     set(gca,'xtick',1:10000:length(timeSplit))
%      set(gca,'xTickLabel',timeLab)
%       set(gca,'FontSize',20);
%      set(gca,'YDir', 'normal');
    set(gca,'xtick',[],'ytick',[])
%   xlabel('Time (s)');
% ylabel('Trials');
end
%% Channel psth
figure;
fireRate = []
for iChan = sigChannelSpike
    spikeRast = spikeRastChan{iChan};
    tSpike = [];
    for iTrial = 1:size(spikeRast,1)
        tSpikeTrial = find(logical(spikeRast(iTrial,:)));
        tSpike = [tSpike tSpikeTrial];
        fireRate(iChan,iTrial,:) = psth(tSpikeTrial,40, fs,size(spikeRast,1),size(spikeRast,2),0);
    end
    subplot(size(chanMap,1),size(chanMap,2),find(ismember(chanMap',iChan)));
    %figure;
   psth(tSpike,20, fs,size(spikeRast,1),size(spikeRast,2),1);
%    set(gca,'FontSize',20);
%    set(gca,'xTickLabel',timeLab)
%     ylim([0 10])
%     xlabel('Time (s)');
% ylabel('Spikes/s');
    set(gca,'xtick',[],'ytick',[])
end
%% Channel plot Spikes
figure;
for iChan = sigChannelSpike
    spikeSplitChan = spikeSplitNeg{iChan};
    timeSpike = linspace(0, 2*0.0025, size(spikeSplitChan,2));
    %figure;
    subplot(size(chanMap,1),size(chanMap,2),find(ismember(chanMap',iChan)));
    plot(timeSpike.*1000,spikeSplitChan,'color',[0 0 0] +0.75);
    hold on;
    plot(timeSpike.*1000,mean(spikeSplitChan,1),'color',[0 0 0]);
     set(gca,'FontSize',20);
%      xlabel('Time (ms)')
%      ylabel('Spike (\muV)');
     xlim([0 5]);
%     title(['Channel ' num2str(iChan)]);
    ylim([-20 20]);
    set(gca,'xtick',[],'ytick',[])
end

%%


spikeTimeId = ones(size(timeSpikeChan));
%spikeTimeId(timeSpikeChan>-0.5&timeSpikeChan<=0)=2;
spikeTimeId(timeSpikeChan>0)=2;
%spikeTimeId(timeSpikeChan>0.5&timeSpikeChan<=1)=4;
[coeff,score,latent,tsquared,explained,mu] = pca(spikeSplitChan);
figure;
scatter(score(:,1),score(:,2),[],spikeTimeId,'filled');
colormap(jet(4096));


figure;
scatter3(score(:,1),score(:,2),score(:,3),[],spikeTimeId,'filled');
colormap(jet(4096));

figure;
plot(coeff(:,1));
hold on;
plot(coeff(:,2));
hold on;
plot(coeff(:,3));
